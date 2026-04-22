"""Python package entrypoint for mathcore_for_cv bindings."""

from __future__ import annotations

HAS_WRAPPER = False
WRAPPER_IMPORT_ERROR: ImportError | None = None

try:
    from .mathcore_for_cv import *  # noqa: F401,F403
except ImportError as exc:
    WRAPPER_IMPORT_ERROR = exc
else:
    HAS_WRAPPER = True
