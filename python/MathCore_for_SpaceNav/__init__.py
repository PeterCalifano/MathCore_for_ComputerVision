"""Python package entrypoint for MathCore_for_SpaceNav bindings."""

from __future__ import annotations

HAS_WRAPPER = False
WRAPPER_IMPORT_ERROR: ImportError | None = None

try:
    from .MathCore_for_SpaceNav import *  # noqa: F401,F403
except ImportError as exc:
    WRAPPER_IMPORT_ERROR = exc
else:
    HAS_WRAPPER = True
