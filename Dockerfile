FROM mambaorg/micromamba:latest

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    hmmer=3.4 \
    python=3.13 \
    pip=26.0.1 \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 2. Set the working directory inside the container
WORKDIR /app

# 3. Copy your requirements and install them
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copy your actual script/logic into the image
COPY needle ./needle
COPY scripts ./scripts

# (Optional) Add your app directory to the PATH
ENV PATH="/app:${PATH}"
ENV PYTHONPATH="/app"

# We don't use ENTRYPOINT or CMD because Nextflow 
# will override them to run its own wrapper script.
