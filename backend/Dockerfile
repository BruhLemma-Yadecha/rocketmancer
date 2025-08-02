FROM python:3.11-slim AS base

# 1. Pin to a known uv image for reproducibility
COPY --from=ghcr.io/astral-sh/uv:python3.11-bookworm-slim /usr/local/bin/uv /usr/local/bin/uvx /usr/local/bin/

ENV PYTHONDONOTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# 2. Copy in just the pyproject manifests plus your rocketmancer package
#    so that 'uv sync' sees the local path dependency
COPY backend/pyproject.toml backend/uv.lock ./
COPY packages/ ../packages/
COPY LICENSE ../LICENSE
COPY README.md ../README.md

# 3. Install all dependencies, including rocketmancer, exactly per uv.lock
RUN uv sync --locked

# 4. Copy the rest of your Django app
COPY backend/ .

EXPOSE 8000

# 5. Launch via uv-wrapped gunicorn
CMD ["uv", "run", "gunicorn", "--bind", "0.0.0.0:8000", "rocketmancer_backend.wsgi:application"]