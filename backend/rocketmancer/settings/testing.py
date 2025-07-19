"""
Testing settings for rocketmancer project.
"""

from .base import *

# Use a simple secret key for testing
SECRET_KEY = "test-secret-key-not-for-production"

# Debug should be False for testing to catch template errors
DEBUG = False

ALLOWED_HOSTS = ["testserver", "localhost"]

# Use in-memory SQLite for faster tests
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": ":memory:",
    }
}

# Disable CORS for testing
CORS_ALLOWED_ORIGINS = []

# Disable logging during tests
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "null": {
            "class": "logging.NullHandler",
        },
    },
    "root": {
        "handlers": ["null"],
    },
}

# Speed up password hashing for tests
PASSWORD_HASHERS = [
    "django.contrib.auth.hashers.MD5PasswordHasher",
]
