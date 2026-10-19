FROM python:3.12-slim

WORKDIR /app

COPY pyproject.toml ./
COPY app ./app
COPY static ./static
RUN pip install --no-cache-dir .

# OZONE build listens on 8002 (legacy 8000, corrected 8001), so the same
# number is correct everywhere -- inside the container, published on the
# host, and in any hosting panel's "internal"/"external" port fields.
# This removes the outside/inside mismatch that otherwise causes a 502
# when a panel assumes the two match.
EXPOSE 8002

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8002"]
