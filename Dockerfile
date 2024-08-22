FROM python:3.12.5
LABEL maintainer=yy5@sanger.ac.uk
LABEL org.opencontainers.image.licenses="MIT"

# Create the working directory
WORKDIR /myworkdir

# Copy source code
COPY src /myworkdir/src/
COPY README.md /myworkdir/src/
COPY pyproject.toml /myworkdir/

# build
RUN pip install --editable .