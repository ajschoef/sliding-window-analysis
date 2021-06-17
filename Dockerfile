FROM python:3.9-slim-buster

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
ENV PYTHONUNBUFFERED=1

# Install dependencies:
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy files
COPY sliding_window sliding_window
COPY data/raw data/raw
COPY run.py .

# Run the application
ENTRYPOINT ["python", "run.py"]
