# Bundles the two compiled analysis tools — freddie-core and putagene — on top of
# a base image that already carries the conda dependencies. Snakemake and tksm2
# are deliberately NOT here: Snakemake orchestrates from the host, and tksm2 /
# simulation is a separate workflow. BASE_IMAGE defaults to the Seqera Containers
# image built from docker/seqera.yaml:
#
#   docker build -f docker/Dockerfile -t scfreddie:latest .
#
# Regenerate the base when seqera.yaml changes (https://seqera.io/containers) and
# update the BASE_IMAGE tag below.
#
# No Seqera image? Build the env yourself from a plain micromamba base:
#   docker build -f docker/Dockerfile \
#     --build-arg BASE_IMAGE=mambaorg/micromamba:1.5.8-jammy --build-arg SETUP_ENV=true \
#     -t scfreddie:latest .
#
# HiGHS-only — no Gurobi (the conda package lacks the C header, and the license is
# a hassle in a container). Build with Gurobi locally instead; see the README.
#
#   docker run --rm scfreddie:latest freddie-core --help
#   docker run --rm scfreddie:latest putagene --help

# Seqera Containers image built from docker/seqera.yaml (regenerate + update on dep changes).
ARG BASE_IMAGE=community.wave.seqera.io/library/scfreddie:a1478580f75ad774
FROM ${BASE_IMAGE}

USER root
# SETUP_ENV=true installs environment.yaml into the base env — needed only when the
# base image doesn't already provide it (e.g. plain micromamba, not a Seqera image).
ARG SETUP_ENV=false
ENV SCFREDDIE_HOME=/opt/scfreddie \
    CONDA_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH

COPY environment.yaml /tmp/environment.yaml
RUN if [ "$SETUP_ENV" = "true" ]; then \
        micromamba install -y -n base -f /tmp/environment.yaml && micromamba clean -a -y ; \
    fi

# pigz decompresses the reads in the flexiplex/demux rules. It's in seqera.yaml
# now, but install it here too so an older base still works after a rebuild.
RUN command -v pigz >/dev/null 2>&1 || \
    { micromamba install -y -n base -c conda-forge pigz && micromamba clean -a -y ; }

# scikit-image: scanpy's sc.pp.scrublet needs it for automatic thresholding, but
# scanpy doesn't pull it in. In seqera.yaml now; install here for older bases.
RUN python -c "import skimage" >/dev/null 2>&1 || \
    { micromamba install -y -n base -c conda-forge scikit-image && micromamba clean -a -y ; }

# python-igraph: sc.tl.leiden(flavor="igraph") needs it; scanpy doesn't pull it in.
# In seqera.yaml now; install here for older bases.
RUN python -c "import igraph" >/dev/null 2>&1 || \
    { micromamba install -y -n base -c conda-forge python-igraph && micromamba clean -a -y ; }

# Build just the two analysis binaries via their own CMake (not the superbuild, so
# no tksm2). -j4 keeps peak memory sane; bare -j would spawn unbounded jobs and OOM.
WORKDIR ${SCFREDDIE_HOME}
COPY freddie-core ${SCFREDDIE_HOME}/freddie-core
COPY sclr-putagene ${SCFREDDIE_HOME}/sclr-putagene
RUN cmake -S freddie-core -B freddie-core/build && cmake --build freddie-core/build -j4 && \
    cmake -S sclr-putagene -B sclr-putagene/build && cmake --build sclr-putagene/build -j4 && \
    ln -sf ${SCFREDDIE_HOME}/freddie-core/build/freddie-core /opt/conda/bin/freddie-core && \
    ln -sf ${SCFREDDIE_HOME}/sclr-putagene/build/putagene    /opt/conda/bin/putagene

CMD ["freddie-core", "--help"]
