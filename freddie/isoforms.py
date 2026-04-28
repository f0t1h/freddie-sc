import enum
import functools
from dataclasses import dataclass, field
import os
import pickle

from freddie.ilp import ALN_T_MAP, FredILP, IlpParams, UnsolvableILP, TimeoutILP
from freddie.segment import CanonIntervals, PairedInterval, aln_t
from freddie.split import Interval, Read, Tint

import numpy as np
import pulp


class timeoutStrat(enum.IntEnum):
    stop = 0
    subsample = 1


@dataclass
class IsoformsParams:
    max_isoform_count: int = 20
    min_read_support: int = 3
    timeout_stategy: timeoutStrat = timeoutStrat.stop
    ilp_params: IlpParams = field(default_factory=IlpParams)

    def __post_init__(self):
        assert 1 <= self.max_isoform_count
        assert 1 <= self.min_read_support


@dataclass
class IntervalSupport(Interval):
    support: float = 0.0


@dataclass
class ProblemSize:
    contig: str
    tint_id: int
    read_count: int
    canon_interval_count: int
    matrix_rows: int
    matrix_cols: int
    unique_read_patterns: int
    cell_type_count: int
    row_cell_type_count: int
    intronic_cell_count: int
    intron_run_count: int
    unique_intron_run_count: int
    estimated_binary_vars: int
    estimated_integer_vars: int
    estimated_constraints: int

    @staticmethod
    def header() -> str:
        return "\t".join(ProblemSize.__dataclass_fields__)

    def __str__(self) -> str:
        return "\t".join(str(getattr(self, field)) for field in self.__dataclass_fields__)


@functools.total_ordering
class Isoform:
    """
    Isoform class.

    Attributes:
        tid: Tint ID
        reads: List of reads comprising the isoform
        iid: Isoform index
        contig: Contig
        strand: Strand of the isoform if it can be determined from the read polyA tails (i.e. + or -). Otherwise, "."
        exons: List of genomic intervals (i.e. exons) comprising the isoform with support values.
                The support value is computed by adding the number of bases covered by
                each read in the interval and dividing by the interval length.
        cell_types: Tuple of cell types

    Methods:
        __eq__: Equality operator
        __lt__: Less than operator
        __repr__: String representation of the isoform in GTF format
    """

    def __init__(
        self,
        tid: int,
        contig: str,
        reads: list[Read],
        isoform_index: int,
    ) -> None:
        self.tid = tid
        self.reads = reads
        self.iid = isoform_index
        self.contig = contig
        self.exons: list[IntervalSupport] = list()
        cell_types_set = set()
        for read in self.reads:
            if len(cell_types_set) == 0:
                cell_types_set.add("NA")
            for ct in read.cell_types:
                cell_types_set.add(ct)
        self.cell_types = tuple(sorted(cell_types_set))

        canon_ints = CanonIntervals(self.reads)
        for i in range(10):
            canon_ints.pop(i)
        intervals: list[IntervalSupport] = list()
        for i in canon_ints.intervals:
            if (e_cnt := len(i.exonic_ridxs())) > len(i.intronic_ridxs()):
                intervals.append(IntervalSupport(i.start, i.end, e_cnt))
        intervals.sort()
        for i in intervals:
            # Add first exon
            if len(self.exons) == 0:
                self.exons.append(i)
                continue
            # Current interval is not adjacent to previous interval: add new exon
            if self.exons[-1].end < i.start:
                self.exons.append(i)
                continue
            # Current interval is adjacent to previous interval: merge exons and update support
            e = self.exons[-1]

            self.exons[-1] = IntervalSupport(
                e.start,
                i.end,
                (len(e) * e.support + len(i) * i.support) / (len(e) + len(i)),
            )

        self.strand = "."
        for read in self.reads:
            if read.polyAs[0].length > 0:
                self.strand = "-"
                break
            if read.polyAs[1].length > 0:
                self.strand = "+"
                break

    def __eq__(self, __value: object) -> bool:
        if not isinstance(__value, Isoform):
            return NotImplemented
        return self.contig == __value.contig and self.exons == __value.exons

    def __lt__(self, __value: object) -> bool:
        if not isinstance(__value, Isoform):
            return NotImplemented
        if self.contig != __value.contig:
            # If contig is a number, sort by number
            if self.contig.isnumeric() and __value.contig.isnumeric():
                return int(self.contig) < int(__value.contig)
            return self.contig < __value.contig
        return self.exons < __value.exons

    def __repr__(self) -> str:
        gtf_records = list()
        gene_id = f"{self.contig}_{self.tid}"
        isoform_id = f"{gene_id}_{self.iid}"
        gtf_records.append(
            "\t".join(
                [
                    self.contig,
                    "freddie",
                    "transcript",
                    f"{self.exons[0].start + 1}",
                    f"{self.exons[-1].end}",
                    ".",
                    self.strand,
                    ".",
                    " ".join(
                        [
                            f'gene_id "{gene_id}";',
                            f'transcript_id "{isoform_id}";',
                            f'read_support "{len(self.reads)}";',
                            f'cell_types "{",".join(self.cell_types)}";',
                        ]
                    ),
                ]
            )
        )
        for idx, exon in enumerate(self.exons, start=1):
            gtf_records.append(
                "\t".join(
                    [
                        self.contig,
                        "freddie",
                        "exon",
                        f"{exon.start + 1}",
                        f"{exon.end}",
                        ".",
                        self.strand,
                        ".",
                        " ".join(
                            [
                                f'gene_id "{gene_id}";',
                                f'transcript_id "{isoform_id}";',
                                f'read_support "{exon.support:.2f}";',
                                f'exon_number "{idx}";',
                            ]
                        ),
                    ]
                )
            )
        return "\n".join(gtf_records)


def get_isoforms(
    tint: Tint,
    params: IsoformsParams = IsoformsParams(),
) -> tuple[Tint, list[Isoform]]:
    """
    Returns isoforms for the given Tint.

    Args:
        tint: Tint
        params: Isoform params dataclass

    Returns:
        isoforms: list[Isoform]
    """
    assert params.min_read_support > 0
    reads: list[Read] = tint.reads
    isoforms = list()
    for isoform_index in range(params.max_isoform_count):
        try:
            canon_ints, recycling_bin, isoform_bin, unsampled_bin = run_ilp_loop(
                reads, params
            )
        except (UnsolvableILP, TimeoutILP) as e:
            os.makedirs("tints", exist_ok=True)
            fname = f"tints/{str(e).replace(' ', '')}.contig_{tint.contig}.tint_{tint.tid}.pickle"
            with open(fname, "wb+") as f:
                pickle.dump(reads, f)
            break
        if len(isoform_bin) < params.min_read_support:
            break
        isoform_reads = [reads[ridx] for ridx in isoform_bin]
        recycling_reads = [reads[ridx] for ridx in recycling_bin]
        if len(unsampled_bin) > 0:
            unsampled_isoform_reads, unsampled_recycing_reads = get_compatible_reads_bins(
                [reads[ridx] for ridx in unsampled_bin],
                isoform_reads,
                canon_ints,
                [aln_t.exon] * (len(canon_ints.intervals) - 1),
                params.ilp_params.max_correction_len,
            )
            isoform_reads.extend(unsampled_isoform_reads)
            recycling_reads.extend(unsampled_recycing_reads)

        isoform = Isoform(
            tid=tint.tid,
            contig=tint.contig,
            reads=isoform_reads,
            isoform_index=isoform_index,
        )
        isoforms.append(isoform)
        # Remove the reads that were used to construct the isoform
        reads = recycling_reads
        if len(reads) < params.min_read_support:
            break
    return tint, isoforms


def run_ilp(
    canon_ints: CanonIntervals, params: IsoformsParams
) -> tuple[int, list[int], list[int]]:
    """
    Args:
        canon_ints: CanonIntervals
        params: Isoforms params dataclass

    Returns:
        recycling_bin: List of ridxs belonging to the recycling bin
        isoform_bin: List of ridxs belonging to the isoform bin
    """
    ilp: FredILP = FredILP(canon_ints, params.ilp_params)
    ilp.build_model(K=2)
    status, (recycling_bin, isoform_bin) = ilp.solve()
    return status, recycling_bin, isoform_bin


def canonize_reads(reads):
    canon_ints = CanonIntervals(reads)
    for i in range(10):
        canon_ints.pop(i)
    return canon_ints


def get_problem_size(
    tint: Tint,
    params: IsoformsParams = IsoformsParams(),
) -> ProblemSize:
    """
    Estimate the first ILP size for a Tint without building or solving the model.
    """
    canon_ints = canonize_reads(tint.reads)
    rows, cell_types = _get_ilp_rows_and_cell_types(canon_ints, params.ilp_params)
    matrix = canon_ints.get_matrix()
    K = 2
    M = len(canon_ints.intervals) + 2
    N = len(rows)
    cell_type_set = {cell_type for cts in cell_types for cell_type in cts}
    J = len(cell_type_set)
    row_cell_type_count = sum(len(cts) for cts in cell_types)
    intronic_cell_count = sum(row.count(aln_t.intron) for row in rows)
    intron_runs = [run for row in rows for run in _get_introns(row)]
    intron_run_count = len(intron_runs)
    unique_intron_run_count = len(set(intron_runs))

    binary_vars = _estimate_binary_vars(
        K=K,
        M=M,
        N=N,
        J=J,
        row_cell_type_count=row_cell_type_count,
        intronic_cell_count=intronic_cell_count,
    )
    integer_vars = (K - 1) * (unique_intron_run_count + 1)
    constraints = _estimate_constraints(
        K=K,
        M=M,
        N=N,
        J=J,
        row_cell_type_count=row_cell_type_count,
        intronic_cell_count=intronic_cell_count,
        intron_run_count=intron_run_count,
        unique_intron_run_count=unique_intron_run_count,
    )

    return ProblemSize(
        contig=tint.contig,
        tint_id=tint.tid,
        read_count=len(tint.reads),
        canon_interval_count=len(canon_ints.intervals),
        matrix_rows=matrix.shape[0],
        matrix_cols=matrix.shape[1],
        unique_read_patterns=N,
        cell_type_count=J,
        row_cell_type_count=row_cell_type_count,
        intronic_cell_count=intronic_cell_count,
        intron_run_count=intron_run_count,
        unique_intron_run_count=unique_intron_run_count,
        estimated_binary_vars=binary_vars,
        estimated_integer_vars=integer_vars,
        estimated_constraints=constraints,
    )


def _get_ilp_rows_and_cell_types(
    canon_ints: CanonIntervals,
    params: IlpParams,
) -> tuple[tuple[tuple[aln_t, ...], ...], tuple[tuple[str, ...], ...]]:
    data: dict[tuple[tuple[aln_t, ...], tuple[str, ...]], None] = dict()
    for idx, row in enumerate(canon_ints.get_matrix()):
        if params.ignore_celltype:
            cell_types = ("",)
        else:
            cell_types = canon_ints.reads[idx].cell_types
        first = len(row) - 1
        last = 0
        for j, aln_type in enumerate(row):
            if aln_type in [aln_t.exon, aln_t.polyA]:
                first = min(first, j)
                last = max(last, j)
        assert first <= last
        key = (
            tuple(aln_t.unaln for _ in range(first))
            + tuple(ALN_T_MAP[i] for i in row[first : last + 1])
            + tuple(aln_t.unaln for _ in range(last + 1, len(row)))
        ), cell_types
        data[key] = None
    keys = tuple(data.keys())
    return tuple(r for r, _ in keys), tuple(cts for _, cts in keys)


def _get_introns(row: tuple[aln_t, ...]) -> list[tuple[int, int]]:
    introns: list[tuple[int, int]] = list()
    start: int | None = None
    for idx, aln_type in enumerate(row):
        if aln_type == aln_t.intron and start is None:
            start = idx
        elif aln_type != aln_t.intron and start is not None:
            introns.append((start, idx - 1))
            start = None
    if start is not None:
        introns.append((start, len(row) - 1))
    return introns


def _estimate_binary_vars(
    K: int,
    M: int,
    N: int,
    J: int,
    row_cell_type_count: int,
    intronic_cell_count: int,
) -> int:
    isoform_count = K - 1
    return (
        N * K  # R2I
        + M * isoform_count  # E2I
        + M * N * isoform_count  # E2IR
        + (M - 1) * N * isoform_count  # EXON_CONTIG2IR
        + (M - 1) * isoform_count  # EXON_CONTIG2I
        + (M + 2) * isoform_count  # C2I
        + M * N * isoform_count  # C2IR
        + (M + 1) * isoform_count  # CHANGE2I
        + J * isoform_count  # I2T
        + M * row_cell_type_count * isoform_count  # C2IRT
        + M * J * isoform_count  # C2IT
        + intronic_cell_count * isoform_count  # OBJ
    )


def _estimate_constraints(
    K: int,
    M: int,
    N: int,
    J: int,
    row_cell_type_count: int,
    intronic_cell_count: int,
    intron_run_count: int,
    unique_intron_run_count: int,
) -> int:
    isoform_count = K - 1
    return (
        N  # read assignment
        + M * N * isoform_count  # E2IR definitions
        + M * isoform_count * (N + 1)  # E2I max
        + (M - 1) * N * isoform_count  # EXON_CONTIG2IR definitions
        + (M - 1) * isoform_count * (N + 1)  # EXON_CONTIG2I max
        + 2 * isoform_count  # C2I boundary definitions
        + M * N * isoform_count  # C2IR definitions
        + M * isoform_count * (N + 1)  # C2I max
        + 4 * (M + 1) * isoform_count  # CHANGE2I xor
        + isoform_count * (row_cell_type_count + J)  # I2T max
        + M * row_cell_type_count * isoform_count  # C2IRT definitions
        + M * isoform_count * (row_cell_type_count + J)  # C2IT max
        + unique_intron_run_count * isoform_count  # GAPI definitions
        + isoform_count  # CHANGE2I_sum definition
        + 3 * (M - 1) * isoform_count  # exon contiguity and
        + 3 * M * J * isoform_count  # cell type and
        + isoform_count  # polyA constraint
        + intron_run_count * isoform_count  # gap constraints
        + 3 * intronic_cell_count * isoform_count  # OBJ and
    )


def run_ilp_loop(
    reads: list[Read],
    params: IsoformsParams,
) -> tuple[CanonIntervals, list[int], list[int], list[int]]:
    """
    Run ILP with the given reads and return the recycling reads and isoform reads.
    If the ILP fails to find an optimal solution, iteratively keep halving the
    number of reads until an optimal solution is found or the number of reads
    drops below the minimum read support.

    Args:
        reads: List of reads
        params: IsoformsParams object

    Returns:
        canon_ints: Canonical intervals of the recycling + isoform reads
        recycling_ridxs: List of recycling read indices
        isoform_ridxs: List of isoform read indices
        unsampled_ridxs: List of unsampled read indices
    """
    N = len(reads)
    reads_idxs = list(range(N))
    canon_ints = canonize_reads(reads)
    if len(reads) < params.min_read_support:
        return canon_ints, reads_idxs, list(), list()
    try:
        status, recycling_bin, isoform_bin = run_ilp(canon_ints, params)
        assert status == pulp.LpStatusOptimal
        return canon_ints, recycling_bin, isoform_bin, list()
    except TimeoutILP as e:
        timeout_exc = e

    if params.timeout_stategy == timeoutStrat.stop:
        raise timeout_exc
    elif params.timeout_stategy == timeoutStrat.subsample:
        N2 = min(
            int(params.ilp_params.timeLimit * 30),  # Each 30 reads take ~1sec
            N // 2,
        )
    else:
        raise ValueError(f"Invalid timeout strategy: {params.timeout_stategy}")

    sample_ridxs: list[int] = list()
    unsampled_ridxs: list[int] = list()
    S: set[int] = set(np.random.choice(reads_idxs, size=N2, replace=False))
    for ridx in reads_idxs:
        if ridx in S:
            sample_ridxs.append(ridx)
        else:
            unsampled_ridxs.append(ridx)

    _, sub_recycling_bin, sub_isoform_bin, sub_unsampled_bin = run_ilp_loop(
        [reads[ridx] for ridx in sample_ridxs],
        params,
    )
    recycling_bin = [sample_ridxs[ridx] for ridx in sub_recycling_bin]
    isoform_bin = [sample_ridxs[ridx] for ridx in sub_isoform_bin]
    unsampled_ridxs.extend([sample_ridxs[ridx] for ridx in sub_unsampled_bin])

    return canon_ints, recycling_bin, isoform_bin, unsampled_ridxs


def get_compatible_reads_bins(
    unsampled_reads: list[Read],
    isoform_reads: list[Read],
    canon_ints: CanonIntervals,
    i_vect: list[aln_t],
    slack: int,
):
    """
    Split unsampled reads into in/compatible lists of reads. The method is used when the ILP
    was run with less than all reads. A read is compatible if it shares an exon with the isoform
    and it does not add any exons to the isoform. Additionally, the number

    Args:
        unsampled_reads: List of unsampled reads
        isoform_reads: List of isoform reads
        intervals: Canonical intervals
        isoform_vect: Isoform vector
        slack: Slack

    Returns:
        incompatible_reads: List of incompatible reads
        compatible_reads: List of compatible reads
    """
    isoform_intervals = [
        PairedInterval(
            target=Interval(canon_ints.intervals[j].start, canon_ints.intervals[j].end)
        )
        for j, e in enumerate(i_vect[1:-1])
        if e == aln_t.exon
    ]
    isoform_read = Read(
        idx=-1,
        name="",
        strand="",
        intervals=isoform_intervals,
        qlen=0,
        polyAs=(
            Read.PolyA(overhang=0, length=i_vect[0] == aln_t.exon, slack=0),
            Read.PolyA(overhang=0, length=i_vect[-1] == aln_t.exon, slack=0),
        ),
        cell_types=tuple({ct for read in isoform_reads for ct in read.cell_types}),
    )

    compatible_reads: list[Read] = list()
    incompatible_reads: list[Read] = list()
    for read in unsampled_reads:
        cints = CanonIntervals([isoform_read, read])
        for i in range(10):
            cints.pop(i)
        M = cints.get_matrix()
        i_row = M[0, :]
        r_row = M[1, :]
        is_compat = True
        for j in range(M.shape[1]):
            if i_row[j] == aln_t.exon and r_row[j] != aln_t.exon:
                is_compat = False
                break
            if r_row[j] == aln_t.exon and i_row[j] != aln_t.intron:
                L = cints.intervals[j - 1].end - cints.intervals[j - 1].start
                if L > slack:
                    is_compat = False
                    break
        if is_compat:
            compatible_reads.append(read)
        else:
            incompatible_reads.append(read)
    return incompatible_reads, compatible_reads
