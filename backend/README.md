# Backend Development

## Setup
```bash
uv sync
```

## Development Commands
```bash
# Run server
uv run python manage.py runserver

# Tests
uv run pytest
uv run pytest --cov=. --cov-report=html

# Code quality
uv run ruff check .
uv run ruff check --fix .
uv run ruff format .
uv run mypy .

# Django commands
uv run python manage.py migrate
uv run python manage.py makemigrations
uv run python manage.py shell
```

## Make Commands
Use `make help` to see all available commands.