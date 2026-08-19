"""ASGI entrypoint under a filename deployment platforms scan for.

Vercel looks for `app.py`/`index.py`/`server.py`/`main.py`/`wsgi.py`/`asgi.py`, so this
re-exports the application instead of configuring a path to `app/api.py`. It builds
nothing: the app, its routes and its middleware are constructed once, in app/api.py.

    uvicorn asgi:app          # equivalent to uvicorn app.api:app
"""

from __future__ import annotations

from app.api import app

__all__ = ["app"]
