FROM python:3.12-slim

WORKDIR /srv
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app app
COPY data/organism.json data/organism.json

# One worker: the CPU-second rate limiter and the codon tables are per-process,
# and determinism is asserted per process.
CMD ["uvicorn", "app.api:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]
