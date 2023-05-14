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

USER mambauser

# Make directories

RUN mkdir -p /home/mambauser/Code

WORKDIR /home/mambauser/Code

# Pull in the latest version of the code
# Layer RMG-Py
RUN git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
# Layer RMG-database
RUN git clone https://github.com/ReactionMechanismGenerator/RMG-database.git

WORKDIR /home/mambauser/Code/RMG-Py

# Install RMG-Py
RUN micromamba create -y -f environment.yml && \
    micromamba install -n rmg_env -c conda-forge ncurses && \
    micromamba clean --all -f -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=rmg_env

# Activate the environment
ENV PATH /opt/conda/envs/rmg_env/bin:$PATH
ENV PYTHONPATH /home/mambauser/Code/RMG-Py:$PYTHONPATH
ENV PATH /home/mambauser/Code/RMG-Py:$PATH
ENV LD_LIBRARY_PATH /opt/conda/envs/rmg_env/lib
ENV PYTHON=/opt/conda/envs/rmg_env/bin/python
ENV PYCALL_JL_RUNTIME_PYTHON=/opt/conda/envs/rmg_env/bin/python


RUN echo "export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> "$BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh" \
    && echo "export LD_LIBRARY_PATH=/opt/conda/envs/rmg_env/lib" >> "$BASE/envs/rmg_env/etc/conda/activate.d/env_vars.sh" \
    && echo "export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}" >> "$BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh" \
    && echo "unset OLD_LD_LIBRARY_PATH" >> "$BASE/envs/rmg_env/etc/conda/deactivate.d/env_vars.sh" \
    && echo "export PYTHONPATH=/home/mambauser/Code/RMG-Py" >> ~/.bashrc \
    && echo "export PATH=/home/mambauser/Code/RMG-Py:$PATH" >> ~/.bashrc \
    && make \
    && ln -sf /opt/conda/envs/rmg_env/bin/python /opt/conda/envs/rmg_env/bin/python-jl 
RUN  /opt/conda/envs/rmg_env/share/julia/site/julia-1.8.5/bin/julia -e 'using Pkg; \
              Pkg.add(PackageSpec(name="ReactionMechanismSimulator", rev="main")); \
              using ReactionMechanismSimulator'\
    && /opt/conda/envs/rmg_env/bin/python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"


# # Add alias to bashrc - rmge to activate the environment
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

WORKDIR /home/mambauser/Code

RUN git clone https://github.com/ReactionMechanismGenerator/ARC.git

WORKDIR /home/mambauser/Code/ARC

ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/ARC"
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/AutoTST"
ENV PYTHONPATH="${PYTHONPATH}:/home/mambauser/Code/TS-GCN"
ENV PATH /home/mambauser/Code/ARC:$PATH

RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y && \
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


SHELL ["/bin/bash", "-c"]
