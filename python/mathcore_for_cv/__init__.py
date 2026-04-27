"""Python package entrypoint for mathcore_for_cv bindings."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

HAS_WRAPPER = False
WRAPPER_IMPORT_ERROR: Exception | None = None


def _export_public_symbols(module: object) -> None:
    public_names = getattr(module, "__all__", None)
    if public_names is None:
        public_names = [name for name in dir(module) if not name.startswith("_")]

    globals().update({name: getattr(module, name) for name in public_names})


def _load_wrapper_from_build_directory() -> bool:
    try:
        from ._wrapper_build import WRAPPER_LIBRARY_DIRS, WRAPPER_MODULE_PATH
    except ImportError:
        return False

    module_path = Path(WRAPPER_MODULE_PATH)
    if not module_path.exists():
        return False

    for library_dir in WRAPPER_LIBRARY_DIRS:
        if library_dir and library_dir not in sys.path:
            sys.path.insert(0, library_dir)

    spec = importlib.util.spec_from_file_location(f"{__name__}.mathcore_for_cv", module_path)
    if spec is None or spec.loader is None:
        return False

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module

    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(spec.name, None)
        raise

    _export_public_symbols(module)
    return True

try:
    from .mathcore_for_cv import *  # noqa: F401,F403
except ImportError as exc:
    WRAPPER_IMPORT_ERROR = exc
    try:
        HAS_WRAPPER = _load_wrapper_from_build_directory()
    except Exception as build_exc:
        WRAPPER_IMPORT_ERROR = build_exc
    else:
        if HAS_WRAPPER:
            WRAPPER_IMPORT_ERROR = None
else:
    HAS_WRAPPER = True
