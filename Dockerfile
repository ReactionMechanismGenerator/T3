# Base image: ARC (built on mambaorg/micromamba, Ubuntu 24.04).
# Contract inherited from the base image:
#   - non-root user `mambauser`, home `/home/mambauser`, code under `/home/mambauser/Code`
#   - micromamba at MAMBA_ROOT_PREFIX=/opt/conda, MAMBA_DOCKERFILE_ACTIVATE=1
#   - entrypoint `/usr/local/bin/entrywrapper.sh`, which starts as root and drops to mambauser via `runuser`
FROM --platform=linux/amd64 laxzal/arc:latest

# Build the T3 layers as the unprivileged base user
USER mambauser

# Clone the T3 repository into the base image's Code directory
WORKDIR /home/mambauser/Code
RUN git clone --depth 1 --single-branch -b main https://github.com/ReactionMechanismGenerator/T3.git
WORKDIR /home/mambauser/Code/T3

# Create the T3 environment (python 3.14, see environment.yml) and slim the image.
# Cleanups are restricted to paths that always exist so a missing path can't break the && chain.
RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y && \
    rm -rf /home/mambauser/.cache/pip /home/mambauser/.cache/yarn && \
    find /home/mambauser/Code/T3 -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Convenience aliases for interactive shells (not required for the image to run)
RUN echo "export arc_path='/home/mambauser/Code/ARC/'" >> ~/.bashrc && \
    echo "export t3_path='/home/mambauser/Code/T3/'" >> ~/.bashrc && \
    echo "alias t3code='cd \$t3_path'" >> ~/.bashrc && \
    echo "alias t3e='micromamba activate t3_env'" >> ~/.bashrc && \
    echo "alias t3='python /home/mambauser/Code/T3/T3.py input.yml'" >> ~/.bashrc

# Activate the T3 environment for subsequent build steps / interactive use
WORKDIR /home/mambauser/
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=t3_env

# Restore root for the inherited entrypoint (/usr/local/bin/entrywrapper.sh drops to mambauser via runuser)
USER root
