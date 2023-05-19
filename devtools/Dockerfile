# Parent Image
FROM --platform=linux/amd64 mambaorg/micromamba:1.3.1-jammy 

# Set the user to root
USER root
# (otherwise python will not be found)
ARG MAMBA_ROOT_PREFIX=/opt/conda

# Install system dependencies
#
# List of deps and why they are needed:
#  - make, gcc, g++ for building RMG
#  - git for downloading RMG respoitories
#  - wget for downloading conda install script
#  - libxrender1 required by RDKit
RUN apt-get update
RUN apt-get install -y make gcc wget git g++ libxrender1
RUN apt-get clean

# Install conda
ENV PATH="$PATH:/opt/conda/bin"

# Set Bash as the default shell for following commands
SHELL ["/bin/bash", "-c"]

# cd
WORKDIR /home/mambauser
RUN mkdir -p Code
WORKDIR /home/mambauser/Code

# Clone the RMG base and database repositories
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-Py.git
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-database.git

# build the conda environment
WORKDIR /home/mambauser/Code/RMG-Py
RUN micromamba env create --file environment.yml -y

# This runs all subsequent commands inside the rmg_env conda environment
#
# Analogous to just activating the environment, which we can't actually do here
# since that requires running conda init and restarting the shell (not possible
# in a Dockerfile build script)
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=rmg_env

# Set environment variables as directed in the RMG installation instructions
ENV RUNNER_CWD=/home/mambauser/Code
ENV PYTHONPATH="$RUNNER_CWD/RMG-Py:$PYTHONPATH"
ENV PATH="$RUNNER_CWD/RMG-Py:$PATH"
ENV PYTHON=/opt/conda/envs/rmg_env/bin/python
ENV PATH /opt/conda/envs/rmg_env/bin:$PATH

# Build RMG
RUN make

# Install and link Julia dependencies for RMS
RUN python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()" || true
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator' || true 
RUN ln -sfn $(which python-jl) $(which python)

# RMG-Py should now be installed and ready
RUN python-jl -c "from pyrms import rms" || true

# Installing ARC

WORKDIR /home/mambauser/Code

RUN git clone --depth 1  https://github.com/ReactionMechanismGenerator/ARC.git

WORKDIR /home/mambauser/Code/ARC

ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/ARC"
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/AutoTST"
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/TS-GCN"
ENV PATH /home/mambauser/Code/ARC:$PATH

RUN micromamba create -y -f environment.yml && \
    rm -rf /home/mambauser/.cache/yarn \
    rm -rf /home/mambauser/.cache/pip &&\
    rm -rf /home/mambauser/.cache/pip && \
    find -name '*.a' -delete && \
    find -name '*.pyc' -delete && \
    find -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/scipy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/numpy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/pandas -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages -name '*.pyx' -delete && \
    make clean

RUN echo "export PYTHONPATH=/home/mambauser/Code/ARC:$PYTHONPATH" >> ~/.bashrc \
    && echo "export PATH=/home/mambauser/Code/ARC:$PATH" >> ~/.bashrc 

# Install T3
WORKDIR /home/mambauser/Code

RUN git clone https://github.com/ReactionMechanismGenerator/T3.git

WORKDIR /home/mambauser/Code/T3

RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y && \
    rm -rf /home/mambauser/.cache/yarn \
    rm -rf /home/mambauser/.cache/pip \
    && find -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/t3_env/lib/python3.7/site-packages/scipy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/t3_env/lib/python3.7/site-packages/numpy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/t3_env/lib/python3.7/site-packages/pandas -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/t3_env/lib/python3.7/site-packages -name '*.pyx' -delete \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 

# Add alias to bashrc - rmge to activate the environment
# These commands are not necessary for the Docker image to run, but they are useful for the user
RUN echo "alias rmge='micromamba activate rmg_env'" >> ~/.bashrc \
    && echo "alias arce='micromamba activate arc_env'" >> ~/.bashrc \
    && echo "alias t3e='micromamba activate t3_env'" >> ~/.bashrc \
    && echo "alias rmg='python-jl /home/mambauser/Code/RMG/rmg.py input.py'" >> ~/.bashrc \
    && echo "alias deact='micromamba deactivate'" >> ~/.bashrc \
    && echo "export rmgpy_path='/home/mambauser/Code/RMG-Py/'" >> ~/.bashrc \
    && echo "export rmgdb_path='/home/mambauser/Code/RMG-database/'" >> ~/.bashrc \
    && echo "export arc_path='/home/mambauser/Code/ARC/'" >> ~/.bashrc \
    && echo "export autotst_path='/home/mambauser/Code/AutoTST/'" >> ~/.bashrc \
    && echo "export tsgcn_path='/home/mambauser/Code/TS-GCN/'" >> ~/.bashrc \
    && echo "alias t3_path='/home/mambauser/Code/T3/'" >> ~/.bashrc \
    && echo "alias arc='python \$arc_path/ARC.py input.yml > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'" >> ~/.bashrc \
    && echo "alias t3 = 'python \$t3_path/T3.py input.yml > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'" >> ~/.bashrc \
    && echo "alias arcrestart='python \$arc_path/ARC.py restart.yml  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'" >> ~/.bashrc \
    && echo "alias restartarc='python \$arc_path/ARC.py restart.yml  > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)'" >> ~/.bashrc \
    && echo "alias rmgcode='cd \$rmgpy_path'" >> ~/.bashrc \
    && echo "alias rmgdb='cd \$rmgdb_path'" >> ~/.bashrc \
    && echo "alias arccode='cd \$arc_path'" >> ~/.bashrc \
    && echo "t3code='cd \$t3_path'" >> ~/.bashrc \
    && echo "alias conda='micromamba'" >> ~/.bashrc \
    && echo "alias mamba='micromamba'" >> ~/.bashrc

# Set the environment to t3_env
ENV ENV_NAME=t3_env

# Set the user back to mambauser
USER mambauser

SHELL [ "/bin/bash", "--login", "-c" ]
