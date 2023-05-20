# The parent image is the base image that the Dockerfile builds upon.
# The RMG installation instructions suggest Anaconda for installation by source, however, we use micromamba for the Docker image due to its smaller size and less overhead.
FROM --platform=linux/amd64 mambaorg/micromamba:jammy

# Set the user as root
USER root

# Test that the CI runs when this commit is added

# Create a login user named rmguser
ARG NEW_MAMBA_USER=rmguser
ARG NEW_MAMBA_USER_ID=1000
ARG NEW_MAMBA_USER_GID=1000
RUN usermod "--login=${NEW_MAMBA_USER}" "--home=/home/${NEW_MAMBA_USER}" \
        --move-home "-u ${NEW_MAMBA_USER_ID}" "${MAMBA_USER}" && \
    groupmod "--new-name=${NEW_MAMBA_USER}" \
             "-g ${NEW_MAMBA_USER_GID}" "${MAMBA_USER}" && \
    # Update the expected value of MAMBA_USER for the
    # _entrypoint.sh consistency check.
    echo "${NEW_MAMBA_USER}" > "/etc/arg_mamba_user" && \
    :

# Set the environment variables
ARG MAMBA_ROOT_PREFIX=/opt/conda
ENV MAMBA_USER=$NEW_MAMBA_USER
ENV BASE=$MAMBA_ROOT_PREFIX

# Install system dependencies
#
# List of deps and why they are needed:
#  - make, gcc, g++ for building RMG
#  - git for downloading RMG respoitories
#  - wget for downloading conda install script
#  - libxrender1 required by RDKit
# Clean up the apt cache to reduce the size of the image
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    libgomp1\
    libxrender1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Change user to the non-root user
USER $MAMBA_USER

# Make directory for RMG-Py and RMG-database
RUN mkdir -p /home/rmguser/Code

# Change working directory to Code
WORKDIR /home/rmguser/Code

# Clone the RMG base and database repositories. The pulled branches are only the main branches.
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-Py.git
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-database.git

# cd into RMG-Py
WORKDIR /home/rmguser/Code/RMG-Py

# Install RMG-Py and then clean up the micromamba cache
RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y

# Activate the RMG environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=rmg_env

# Set environment variables
# These need to be set in the Dockerfile so that they are available to the build process
ENV PATH /opt/conda/envs/rmg_env/bin:$PATH
ENV PYTHONPATH /home/rmguser/Code/RMG-Py:$PYTHONPATH
ENV PATH /home/rmguser/Code/RMG-Py:$PATH

# Build RMG
RUN make \
    && echo "export PYTHONPATH=/home/rmguser/Code/RMG-Py" >> ~/.bashrc \
    && echo "export PATH=/home/rmguser/Code/RMG-Py:$PATH" >> ~/.bashrc

# Create the link between Julia and RMG-Py. A symbolic link is created between Python-JL and Python.
# We do not use the command in Julia 'using ReactionMechanismSimulator' because the command is only for determining if the package is installed.
# The command does not actually install the package but instead increases the build time, therefore not required.
RUN  python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()" \
    && julia -e 'using Pkg; \
                Pkg.add(PackageSpec(name="ReactionMechanismSimulator", rev="main")); \
                '\
    && ln -sf /opt/conda/envs/rmg_env/bin/python-jl /opt/conda/envs/rmg_env/bin/python

# Compile the RMS during Docker build - This will reduce the time it takes to run RMS for the first time
RUN python -c "from pyrms import rms"

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
    && echo "alias rmg='python-jl /home/rmguser/Code/RMG/rmg.py input.py'" >> ~/.bashrc \
    && echo "alias deact='micromamba deactivate'" >> ~/.bashrc \
    && echo "export rmgpy_path='/home/rmguser/Code/RMG-Py/'" >> ~/.bashrc \
    && echo "export rmgdb_path='/home/rmguser/Code/RMG-database/'" >> ~/.bashrc \
    && echo "export arc_path='/home/rmguser/Code/ARC/'" >> ~/.bashrc \
    && echo "export autotst_path='/home/rmguser/Code/AutoTST/'" >> ~/.bashrc \
    && echo "export tsgcn_path='/home/rmguser/Code/TS-GCN/'" >> ~/.bashrc \
    && echo "alias t3_path='/home/rmguser/Code/T3/'" >> ~/.bashrc \
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

ENV ENV_NAME=t3_env

SHELL [ "/bin/bash", "--login", "-c" ]
