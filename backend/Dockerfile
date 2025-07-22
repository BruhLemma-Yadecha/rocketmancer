FROM python:3.11

ENV PYTHONDONOTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

COPY LICENSE ../LICENSE
COPY README.md ../README.md
COPY backend/pyproject.toml backend/uv.lock ./

RUN uv sync

COPY backend/ .

EXPOSE 8000

CMD ["uv", "run", "gunicorn", "--bind", "0.0.0.0:8000", "rocketmancer.wsgi:application"]
