# Using the a base image of rmgpy:latest from DockerHub - Not the official version but a custom version that is smaller in size
FROM --platform=linux/amd64 laxzal/rmgpy:latest


USER rmguser


# Installing ARC
# Change directory to Code
WORKDIR /home/rmguser/Code

# Clone main branch ARC repository from GitHub and set as working directory
RUN git clone -b main https://github.com/ReactionMechanismGenerator/ARC.git
WORKDIR /home/rmguser/Code/ARC

# Set environment variables for the Docker run and container
ENV PYTHONPATH="${PYTHONPATH}:/home/rmguser/Code/ARC"
ENV PYTHONPATH="${PYTHONPATH}:/home/rmguser/Code/AutoTST"
ENV PYTHONPATH="${PYTHONPATH}:/home/rmguser/Code/TS-GCN"
ENV PATH /home/rmguser/Code/ARC:$PATH

# Install ARC Environment
RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y && \
    rm -rf /home/rmguser/.cache/yarn \
    rm -rf /home/rmguser/.cache/pip &&\
    rm -rf /home/rmguser/.cache/pip && \
    find -name '*.a' -delete && \
    find -name '*.pyc' -delete && \
    find -name '__pycache__' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/scipy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/numpy -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages/pandas -name 'tests' -type d -exec rm -rf '{}' '+' && \
    find /opt/conda/envs/arc_env/lib/python3.7/site-packages -name '*.pyx' -delete && \
    rm -rf /opt/conda/envs/arc_env/lib/python3.7/site-packages/uvloop/loop.c &&\
    make clean

# Install T3

# Change directory to Code
WORKDIR /home/rmguser/Code

# Clone main branch T3 repository from GitHub and set as working directory
RUN git clone -b main https://github.com/ReactionMechanismGenerator/T3.git
WORKDIR /home/rmguser/Code/T3

# Install T3 Environment
RUN micromamba create -y -f environment.yml && \
    micromamba clean --all -f -y && \
    rm -rf /home/rmguser/.cache/yarn \
    rm -rf /home/rmguser/.cache/pip \
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
RUN echo "alias arce='micromamba activate arc_env'" >> ~/.bashrc \
    && echo "alias arc='python /home/rmguser/Code/ARC/ARC.py input.yml'" >> ~/.bashrc \
    && echo "alias deact='micromamba deactivate'" >> ~/.bashrc \
    && echo "export arc_path='/home/rmguser/Code/ARC/'" >> ~/.bashrc \
    && echo "export t3_path='/home/rmguser/Code/T3/'" >> ~/.bashrc \
    && echo "alias arccode='cd \$arc_path'" >> ~/.bashrc \
    && echo "alias t3code='cd \$t3_path'" >> ~/.bashrc \
    && echo "alias t3e='micromamba activate t3_env'" >> ~/.bashrc \
    && echo "alias t3='python /home/rmguser/Code/T3/T3.py input.yml'" >> ~/.bashrc 

SHELL ["/bin/bash", "-c"]
