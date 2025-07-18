.PHONY: help build up down logs clean dev test lint format

# Default target
help:
	@echo "Available commands:"
	@echo "  make build     - Build all Docker containers"
	@echo "  make up        - Start all services"
	@echo "  make down      - Stop all services"
	@echo "  make dev       - Start development environment"
	@echo "  make logs      - Show logs from all services"
	@echo "  make clean     - Clean up containers and volumes"
	@echo "  make test      - Run tests"
	@echo "  make lint      - Run linting"
	@echo "  make format    - Format code"

# Docker commands
build:
	docker-compose build

up:
	docker-compose up -d

down:
	docker-compose down

logs:
	docker-compose logs -f

clean:
	docker-compose down -v --remove-orphans
	docker system prune -f

# Development
dev:
	@echo "Starting development environment..."
	@if [ ! -f .env ]; then \
		echo "Creating .env from .env.example..."; \
		cp .env.example .env; \
	fi
	docker-compose up --build

# Testing
test:
	@echo "Running backend tests..."
	docker-compose exec backend python manage.py test
	@echo "Running frontend tests..."
	docker-compose exec frontend npm test

# Code quality
lint:
	@echo "Linting backend..."
	docker-compose exec backend flake8 .
	@echo "Linting frontend..."
	docker-compose exec frontend npm run lint

format:
	@echo "Formatting backend code..."
	docker-compose exec backend black .
	@echo "Formatting frontend code..."
	docker-compose exec frontend npm run format

# Database
migrate:
	docker-compose exec backend python manage.py migrate

makemigrations:
	docker-compose exec backend python manage.py makemigrations

# Utility
shell-backend:
	docker-compose exec backend python manage.py shell

shell-frontend:
	docker-compose exec frontend sh

install-deps:
	@echo "Installing backend dependencies..."
	docker-compose exec backend pip install -r requirements.txt
	@echo "Installing frontend dependencies..."
	docker-compose exec frontend npm install