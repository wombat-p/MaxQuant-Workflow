FROM nfcore/base
LABEL authors="Veit Schwämmle" \
      description="Docker image containing all software requirements for the nf-core/wombatp pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda install -c conda-forge -n base mamba
RUN mamba env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-wombatp-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-wombatp-1.0dev > nf-core-wombatp-1.0dev.yml

# TEMPORARY SOLUTION FOR SDRF-PIPELINES
RUN conda install pip git
RUN ls
RUN git clone -b extended https://github.com/elixir-proteomics-community/sdrf-pipelines
RUN cd sdrf-pipelines; python setup.py install; cd ..
ENV PYTHON_EGG_CACHE /tmp/
#RUN pip install git+https://github.com/elixir-proteomics-community/sdrf-pipelines.git@extended


# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
