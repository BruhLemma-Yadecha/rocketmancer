#!/bin/bash
# Sync version from VERSION file to all project files

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Get version from parameter or VERSION file
if [ $# -eq 1 ]; then
    VERSION=$1
    echo "$VERSION" > VERSION
else
    # Get version from VERSION file
    if [ ! -f "VERSION" ]; then
        echo -e "${RED}Error: VERSION file not found${NC}"
        exit 1
    fi
    VERSION=$(cat VERSION | tr -d '\n')
fi

echo -e "${YELLOW}Syncing version $VERSION to all project files...${NC}"

# Update frontend package.json
if [ -f "frontend/package.json" ]; then
    echo -e "${GREEN}Updating frontend/package.json...${NC}"
    sed -i "s/\"version\": \"[^\"]*\"/\"version\": \"$VERSION\"/" frontend/package.json
fi

# Update backend pyproject.toml
if [ -f "backend/pyproject.toml" ]; then
    echo -e "${GREEN}Updating backend/pyproject.toml...${NC}"
    sed -i "s/version = \"[^\"]*\"/version = \"$VERSION\"/" backend/pyproject.toml
fi

echo -e "${GREEN}Version sync completed: $VERSION${NC}"
