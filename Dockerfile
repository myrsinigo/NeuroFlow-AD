FROM continuumio/miniconda3:24.1.2-0

WORKDIR /app
COPY environment.yml /app/environment.yml
RUN conda env create -f environment.yml && conda clean -afy
ENV PATH=/opt/conda/envs/ad-sc/bin:$PATH
COPY bin /app/bin
