# Container for the sex-survival pipeline. Every step runs inside it.
#
# Built and pushed once by the maintainer (users don't need this;
# Nextflow pulls the image automatically from Docker Hub):
#   docker buildx build -f envs/Dockerfile --platform linux/amd64,linux/arm64 \
#       -t caglayaniris/sex-survival:0.1.0 --push .
FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/analysis.yml /tmp/analysis.yml

# install the package list, then delete caches to keep the image small
RUN micromamba install -y -n base -f /tmp/analysis.yml && \
    micromamba clean --all --yes

# make the environment's programs visible to Nextflow's commands;
# give matplotlib a writable folder (the container runs as your user)
ENV PATH=/opt/conda/bin:$PATH \
    MPLCONFIGDIR=/tmp/matplotlib
