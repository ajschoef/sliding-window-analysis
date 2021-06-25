FROM python:3.9-slim-buster

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
ENV PYTHONUNBUFFERED=1

# Install dependencies:
COPY requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Copy files
COPY sliding_window sliding_window
ARG path_to_data
COPY $path_to_data $path_to_data
COPY run.py .

# Run the application
ENTRYPOINT ["python", "run.py"]
