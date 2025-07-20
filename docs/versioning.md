# Version Management

This project uses automated semantic versioning with `release-please`.

## Version Files

The single source of truth is the `VERSION` file. All other version references are automatically updated:
- `frontend/package.json`
- `backend/pyproject.toml`

## Manual Sync (if needed)

If version files get out of sync, you can manually sync them:

```bash
./scripts/sync-version.sh
```

Or set a specific version:

```bash
./scripts/sync-version.sh 1.2.3
```
