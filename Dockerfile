FROM python:3.12-slim

WORKDIR /srv
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app app
COPY data/organism.json data/organism.json

# One worker: the CPU-second rate limiter and the codon tables are per-process,
# and determinism is asserted per process.
#
# Shell form on purpose: container hosts inject the port they expect the process to
# listen on (Cloud Run sends PORT=8080 and fails its health check otherwise), and a
# JSON-array CMD would pass "${PORT:-8000}" through as a literal string. `exec` keeps
# uvicorn as PID 1 so it receives SIGTERM and shuts down cleanly on redeploy.
CMD ["sh", "-c", "exec uvicorn app.api:app --host 0.0.0.0 --port ${PORT:-8000} --workers 1"]
