# Getting Started with Development

This guide will help you set up the Rocketmancer Web development environment.

## Prerequisites

- Docker and Docker Compose
- Git
- Node.js 18+ (for local frontend development)
- Python 3.11+ (for local backend development)

## Quick Start

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd rocketmancer-web
   ```

2. **Set up environment variables**
   ```bash
   cp .env.example .env
   # Edit .env with your configuration
   ```

3. **Start the development environment**
   ```bash
   make dev
   ```

4. **Access the application**
   - Frontend: http://localhost
   - Backend API: http://localhost:8000
   - Admin Panel: http://localhost:8000/admin

## Development Workflow

### Using Make Commands

The project includes a Makefile with common development tasks:

```bash
make help          # Show available commands
make build         # Build Docker containers
make up            # Start services
make down          # Stop services
make logs          # View logs
make test          # Run tests
make lint          # Run linting
make clean         # Clean up containers
```

### Local Development (without Docker)

#### Backend Setup

1. Navigate to backend directory:
   ```bash
   cd backend
   ```

2. Create virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements/development.txt
   ```

4. Run migrations:
   ```bash
   python manage.py migrate
   ```

5. Start development server:
   ```bash
   python manage.py runserver
   ```

#### Frontend Setup

1. Navigate to frontend directory:
   ```bash
   cd frontend
   ```

2. Install dependencies:
   ```bash
   npm install
   ```

3. Start development server:
   ```bash
   npm run dev
   ```

## Project Structure

```
rocketmancer-web/
├── backend/                 # Django REST API
│   ├── rocketmancer/       # Django project
│   ├── optimizer/          # Optimization app
│   ├── requirements/       # Python dependencies
│   └── tests/             # Backend tests
├── frontend/               # React frontend
│   ├── src/
│   │   ├── components/    # React components
│   │   ├── services/      # API services
│   │   ├── hooks/         # Custom React hooks
│   │   ├── utils/         # Utility functions
│   │   └── constants/     # App constants
│   └── tests/             # Frontend tests
├── nginx/                  # Nginx configuration
├── docs/                   # Documentation
└── scripts/               # Build/deployment scripts
```

## Testing

### Backend Tests
```bash
# Run all backend tests
make test

# Run specific test file
docker-compose exec backend python manage.py test tests.test_rocket_optimization

# Run with coverage
docker-compose exec backend coverage run --source='.' manage.py test
docker-compose exec backend coverage report
```

### Frontend Tests
```bash
# Run frontend tests
docker-compose exec frontend npm test

# Run with coverage
docker-compose exec frontend npm run test:coverage
```

## Code Quality

### Linting and Formatting

```bash
# Lint all code
make lint

# Format all code
make format

# Backend specific
docker-compose exec backend black .
docker-compose exec backend flake8 .

# Frontend specific
docker-compose exec frontend npm run lint
docker-compose exec frontend npm run format
```

## Troubleshooting

### Common Issues

1. **Port conflicts**: Ensure ports 80 and 8000 are available
2. **Docker issues**: Try `make clean` to reset containers
3. **Permission issues**: Check Docker permissions on your system

### Getting Help

- Check the [FAQ](../FAQ.md)
- Review [Architecture Documentation](../architecture/README.md)
- Open an issue on the project repository