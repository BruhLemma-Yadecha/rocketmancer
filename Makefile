.PHONY: help build up down logs clean dev test lint format lint-backend lint-fix-backend format-backend typecheck-backend test-backend test-backend-cov dev-backend sync-backend

# Default target
help:
	@echo "Available commands:"
	@echo ""
	@echo "Docker commands:"
	@echo "  make build          - Build all Docker containers"
	@echo "  make up             - Start all services"
	@echo "  make down           - Stop all services"
	@echo "  make dev            - Start development environment"
	@echo "  make logs           - Show logs from all services"
	@echo "  make clean          - Clean up containers and volumes"
	@echo ""
	@echo "Development commands:"
	@echo "  make test           - Run all tests"
	@echo "  make lint           - Run linting on all code"
	@echo "  make format         - Format all code"
	@echo ""
	@echo "Backend development (uv):"
	@echo "  make sync-backend      - Sync backend dependencies with uv"
	@echo "  make dev-backend       - Start Django development server"
	@echo "  make test-backend      - Run backend tests with pytest"
	@echo "  make test-backend-cov  - Run backend tests with coverage"
	@echo "  make lint-backend      - Run ruff check on backend"
	@echo "  make lint-fix-backend  - Run ruff check with auto-fix"
	@echo "  make format-backend    - Format backend code with ruff"
	@echo "  make typecheck-backend - Type check backend with mypy"

# Docker commands
build:
	docker compose build

up:
	docker compose up -d

down:
	docker compose down

logs:
	docker compose logs -f

clean:
	docker compose down -v --remove-orphans
	docker system prune -f

# Development
dev:
	@echo "Starting development environment..."
	@if [ ! -f .env ]; then \
		echo "Creating .env from .env.example..."; \
		cp .env.example .env; \
	fi
	docker compose up --build

# Testing
test:
	@echo "Running backend tests..."
	docker compose exec backend python manage.py test
	@echo "Running frontend tests..."
	docker compose exec frontend npm test

# Code quality
lint:
	@echo "Linting backend with ruff..."
	cd backend && uv run ruff check .
	@echo "Type checking backend with mypy..."
	cd backend && uv run mypy .
	@echo "Linting frontend..."
	docker compose exec frontend npm run lint

format:
	@echo "Formatting backend code with ruff..."
	cd backend && uv run ruff format .
	@echo "Formatting frontend code..."
	docker compose exec frontend npm run format

# Backend development commands using uv
lint-backend:
	@echo "Running ruff check..."
	cd backend && uv run ruff check .

lint-fix-backend:
	@echo "Running ruff check with auto-fix..."
	cd backend && uv run ruff check --fix .

format-backend:
	@echo "Formatting backend code with ruff..."
	cd backend && uv run ruff format .

typecheck-backend:
	@echo "Type checking backend with mypy..."
	cd backend && uv run mypy .

test-backend:
	@echo "Running backend tests with pytest..."
	cd backend && uv run pytest

test-backend-cov:
	@echo "Running backend tests with coverage..."
	cd backend && uv run pytest --cov=. --cov-report=html --cov-report=term

dev-backend:
	@echo "Starting Django development server..."
	cd backend && uv run python manage.py runserver

sync-backend:
	@echo "Syncing backend dependencies with uv..."
	cd backend && uv sync

# Database
migrate:
	docker compose exec backend python manage.py migrate

makemigrations:
	docker compose exec backend python manage.py makemigrations

# Backend database commands using uv
migrate-backend:
	@echo "Running Django migrations..."
	cd backend && uv run python manage.py migrate

makemigrations-backend:
	@echo "Creating Django migrations..."
	cd backend && uv run python manage.py makemigrations

shell-backend-uv:
	@echo "Starting Django shell..."
	cd backend && uv run python manage.py shell

collectstatic-backend:
	@echo "Collecting static files..."
	cd backend && uv run python manage.py collectstatic --noinput

# Utility
shell-backend:
	docker compose exec backend python manage.py shell

shell-frontend:
	docker compose exec frontend sh

install-deps:
	@echo "Installing backend dependencies with uv..."
	cd backend && uv sync
	@echo "Installing frontend dependencies..."
	docker compose exec frontend npm install