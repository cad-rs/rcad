"""Python API for the open-source RCAD B-rep kernel."""

from rcad._rcad import (
    BRep,
    BooleanOptions,
    TOLERANCE_ABS,
    brep_proj_cylindrical,
    max_face_tolerance,
    resolved_boolean_fuzzy_tol,
    __version__,
)

__all__ = [
    "BRep",
    "BooleanOptions",
    "TOLERANCE_ABS",
    "brep_proj_cylindrical",
    "max_face_tolerance",
    "resolved_boolean_fuzzy_tol",
    "__version__",
]
