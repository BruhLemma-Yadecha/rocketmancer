# Architecture Overview

This document provides a high-level overview of the Rocketmancer Web system architecture.

## System Architecture

Rocketmancer Web follows a modern web application architecture with clear separation between frontend and backend components.

### Components

- **Frontend**: React-based single-page application
- **Backend**: Django REST API with optimization engine
- **Database**: SQLite
- **Web Server**: Nginx (reverse proxy and static file serving)
- **Containerization**: Docker and Docker Compose

## Architecture Diagrams

- [System Overview](system-overview.md)
- [Data Flow](data-flow.md)
- [Component Interactions](component-interactions.md)

## Design Principles

### Separation of Concerns
- Frontend handles UI/UX and user interactions
- Backend manages business logic and data processing
- Database stores persistent data
- Web server handles routing and static assets

### Scalability
- Stateless API design
- Containerized deployment
- Horizontal scaling capabilities
- Caching strategies

### Security
- Input validation and sanitization
- CORS configuration
- HTTPS enforcement

## Technology Stack

### Frontend
- **React**: UI framework
- **Vite**: Build tool and development server
- **Axios**: HTTP client for API communication
- **React Router**: Client-side routing

### Backend
- **Django**: Web framework
- **Django REST Framework**: API framework
- **NumPy/SciPy**: Scientific computing for optimization
- **SQLite**: Database

### Infrastructure
- **Docker**: Containerization
- **Nginx**: Web server and reverse proxy
- **Docker Compose**: Multi-container orchestration

## Data Architecture

### Database Schema
- User management
- Rocket configuration data
- Optimization results
- System configuration

### API Design
- RESTful endpoints
- JSON data format
- Consistent response structure
- Proper HTTP status codes

For detailed architectural information, see the specific architecture documents linked above.