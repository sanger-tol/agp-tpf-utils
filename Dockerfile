FROM python:3.12
LABEL maintainer=yy5@sanger.ac.uk
LABEL org.opencontainers.image.licenses="MIT"

ARG work_dir="/myworkdir"
WORKDIR $work_dir

# Copy source code
COPY src $work_dir/src/
COPY pyproject.toml $work_dir/

RUN pip install --no-cache-dir --upgrade pip \
  && pip install --no-cache-dir . \
  && rm -r $work_dir/*
COPY README.md $work_dir/
