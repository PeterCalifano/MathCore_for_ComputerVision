"""Public Python package entrypoint for mathcore_for_cv."""

from __future__ import annotations

from importlib import import_module, util
from pathlib import Path
import os
import sys
from types import ModuleType

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


def _export_wrapper_module(module_: ModuleType) -> None:
    public_names_ = getattr(module_, "__all__", None)
    if public_names_ is None:
        public_names_ = [name_ for name_ in dir(module_) if not name_.startswith("_")]
    for name_ in public_names_:
        globals()[name_] = getattr(module_, name_)


def _import_build_linked_wrapper() -> ModuleType:
    try:
        from . import _wrapper_build
    except ImportError as exc:
        raise ImportError("No build-linked wrapper metadata is available.") from exc

    module_path_ = Path(_wrapper_build.WRAPPER_MODULE_PATH)
    if not module_path_.is_file():
        raise ImportError(f"Build-linked wrapper module was not found at '{module_path_}'.")

    if os.name == "nt":
        for dll_dir_ in getattr(_wrapper_build, "WRAPPER_LIBRARY_DIRS", []):
            dll_path_ = Path(dll_dir_)
            if dll_path_.is_dir():
                os.add_dll_directory(str(dll_path_))

    package_name_ = __name__.split(".")[-1]
    module_name_ = f"{__name__}.{package_name_}"
    spec_ = util.spec_from_file_location(module_name_, module_path_)
    if spec_ is None or spec_.loader is None:
        raise ImportError(f"Could not load build-linked wrapper spec from '{module_path_}'.")

    module_ = util.module_from_spec(spec_)
    sys.modules[module_name_] = module_
    spec_.loader.exec_module(module_)
    return module_


package_name_ = __name__.split(".")[-1]
wrapper_module_: ModuleType | None = None

try:
    wrapper_module_ = import_module(f".{package_name_}", __name__)
except ImportError:
    try:
        wrapper_module_ = _import_build_linked_wrapper()
    except ImportError as exc:
        WRAPPER_IMPORT_ERROR = exc

if wrapper_module_ is not None:
    _export_wrapper_module(wrapper_module_)
    HAS_WRAPPER = True
