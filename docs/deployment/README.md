# Deployment Guide

This guide covers deploying Rocketmancer Web to various environments.

## Deployment Options

- [Docker Deployment](docker.md)
- [Cloud Deployment](cloud.md)
- [Production Setup](production.md)
- [Environment Configuration](environment.md)

## Quick Deploy

### Using Docker Compose

```bash
# Start the application
docker-compose up -d

# Or use the Makefile
make up
```

### Environment Variables

Copy and configure the environment file:
```bash
cp .env.example .env
# Edit .env with your configuration
```

## Prerequisites

- Docker and Docker Compose
- Git

## Services

The application runs the following services:
- **Backend**: Django API server on port 8000
- **Frontend**: React app served via Nginx on port 80
- **Database**: SQLite (included with Django)

## Security Considerations

- Use HTTPS in production
- Configure proper firewall rules
- Enable CORS properly
- Use environment variables for secrets
- Secure the Django secret key

## Monitoring

- Monitor application logs via `docker-compose logs`
- Check container health with `docker-compose ps`
- Monitor resource usage with `docker stats`

## Backup Strategy

- SQLite database backup (backend/db.sqlite3)
- Static file backups
- Configuration backups (.env file)

For detailed deployment instructions, see the specific deployment guides linked above.