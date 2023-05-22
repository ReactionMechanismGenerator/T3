FROM --platform=linux/amd64 mambaorg/micromamba:1.3.1-jammy


#Docker is built on layers - Each Run is a Layer

# Install dependencies

USER root
  # (otherwise python will not be found)
ARG MAMBA_ROOT_PREFIX=/opt/conda
ENV PYTHONDONTWRITEBYTECODE=true
ENV BASE=$MAMBA_ROOT_PREFIX

RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    libgomp1\
    libxrender1 \
    libtinfo6 \
    # Clean up
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Change user to mambauser (Non-root) - Defined in the base image.
# The user can be changed via extra arguments to docker build
# https://github.com/mamba-org/micromamba-docker#changing-the-user-id-or-name
USER mambauser

# Make directory and set as working directory
RUN mkdir -p /home/mambauser/Code
WORKDIR /home/mambauser/Code

# Pull in the main branch of the codes
# Layer RMG-Py
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-Py.git
# Layer RMG-database
RUN git clone -b main https://github.com/ReactionMechanismGenerator/RMG-database.git



# Install RMG-Py
WORKDIR /home/mambauser/Code/RMG-Py
RUN micromamba create -y -f environment.yml

# Activate the environment - This is required by the Micromamba base image
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=rmg_env

# Set environment variables for the Docker run and container
ENV PATH /opt/conda/envs/rmg_env/bin:$PATH
ENV PYTHONPATH /home/mambauser/Code/RMG-Py
ENV PATH /home/mambauser/Code/RMG-Py:$PATH
ENV LD_LIBRARY_PATH /opt/conda/envs/rmg_env/lib
ENV PYTHON=/opt/conda/envs/rmg_env/bin/python
ENV PYCALL_JL_RUNTIME_PYTHON=/opt/conda/envs/rmg_env/bin/python

# Compile RMG-Py
RUN make 

# Install RMS
# The extra arguments are required to install PyCall and RMS in this Dockerfile. Will not work without them.
RUN  /opt/conda/envs/rmg_env/bin/julia -e 'using Pkg; \
            ENV["CONDA_JL_HOME"]="/opt/conda/envs/rmg_env/bin"; \
            ENV["CONDA_JL_CONDA_EXE"]="/usr/bin/micromamba"; \
            Pkg.add("PyCall"); \
              Pkg.build("PyCall"); \
              using PyCall; \
              Pkg.add(PackageSpec(name="ReactionMechanismSimulator", rev="main")); \
              using ReactionMechanismSimulator'\
    && /opt/conda/envs/rmg_env/bin/python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()" 

# Compile RMS
RUN python-jl -c "from pyrms import rms"

# Add aliases to bashrc
RUN echo "alias rmge='micromamba activate rmg_env'" >> ~/.bashrc \
    && echo "alias rmg='python-jl /home/mambauser/Code/RMG/rmg.py input.py'" >> ~/.bashrc \
    && echo "alias t3e='micromamba activate t3_env'" >> ~/.bashrc \
    && echo "alias t3='python /home/mambauser/Code/T3/t3.py input.yml'" >> ~/.bashrc \
    && echo "alias arc='python /home/mambauser/Code/ARC/arc.py input.yml'" >> ~/.bashrc \
    && echo "alias arce='micromamba activate arc_env'" >> ~/.bashrc \
    && echo "alias deact='micromamba deactivate'" >> ~/.bashrc \
    && echo "alias jupyter='jupyter notebook'" >> ~/.bashrc \
    && echo "export rmgpy_path='/home/mambauser/Code/RMG-Py/'" >> ~/.bashrc \
    && echo "export rmgdb_path='/home/mambauser/Code/RMG-database/'" >> ~/.bashrc \
    && echo "export arc_path='/home/mambauser/Code/ARC/'" >> ~/.bashrc \
    && echo "export t3_path='/home/mambauser/Code/T3/'" >> ~/.bashrc \
    && echo "alias rmgcode='cd \$rmgpy_path'" >> ~/.bashrc \
    && echo "alias rmgdb='cd \$rmgdb_path'" >> ~/.bashrc \
    && echo "alias arcode='cd \$arc_path'" >> ~/.bashrc \
    && echo "alias t3code='cd \$t3_path'" >> ~/.bashrc \
    && echo "alias conda='micromamba'" >> ~/.bashrc \
    && echo "alias mamba='micromamba'" >> ~/.bashrc 

# Installing ARC

# Change directory to Code
WORKDIR /home/mambauser/Code

# Clone main branch ARC repository from GitHub and set as working directory
RUN git clone -b main https://github.com/ReactionMechanismGenerator/ARC.git
WORKDIR /home/mambauser/Code/ARC

# Set environment variables for the Docker run and container
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/ARC"
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/AutoTST"
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/TS-GCN"
ENV PATH /home/mambauser/Code/ARC:$PATH

# Install ARC Environment
RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y && \
    rm -rf /home/mambauser/.cache/yarn \
    rm -rf /home/mambauser/.cache/pip &&\
    rm -rf /home/mambauser/.cache/pip && \
    find -name '*.a' -delete && \
    find -name '*.pyc' -delete && \
    rm -rf /opt/conda/envs/arc_env/conda-meta && \
    rm -rf /opt/conda/envs/arc_env/include && \
    find -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/scipy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/numpy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/pandas -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages -name '*.pyx' -delete && \
    rm -rf /opt/conda/envs/arc_env/lib/python3.7/site-packages/uvloop/loop.c &&\
    make clean

# Install T3

# Change directory to Code
WORKDIR /home/mambauser/Code

# Clone main branch T3 repository from GitHub and set as working directory
RUN git clone -b main https://github.com/ReactionMechanismGenerator/T3.git
WORKDIR /home/mambauser/Code/T3

# Install T3 Environment
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


SHELL ["/bin/bash", "-c"]
