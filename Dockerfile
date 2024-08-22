FROM python:3.12.5

# Create the working directory
WORKDIR /myworkdir

# Copy source code
COPY src /myworkdir/src/
COPY README.md /myworkdir/src/
COPY pyproject.toml /myworkdir/
COPY tests /myworkdir/
RUN mkdir -p /myworkdir/venv/ && \
    chmod -R 775 /myworkdir/venv/

# Set up the virtual environment and install dependencies

RUN python3 -m venv --prompt asm-utils venv && \
    pip install --editable .

# Activate virtual environment by default
ENV PATH="/myworkdir/venv/bin:$PATH"

# This step is optional and for debugging purposes
RUN python --version && pip --version