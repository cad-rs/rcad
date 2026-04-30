"""One-shot helper: replace common tolerance scientific literals with crate::tolerance names.

Run from repo root: uv run python rcad/libs/rcad-algorithms/tools/normalize_tolerance_literals.py

Re-run is mostly idempotent (literals already replaced won't match).
"""
from __future__ import annotations

from pathlib import Path

CRATE_ROOT = Path(__file__).resolve().parents[1]
SRC = CRATE_ROOT / "src"
TESTS = CRATE_ROOT / "tests"

# Longest / most specific first (plain substring replace).
#
# Avoid Rust typed floats like `1e-6_f64`: naive replace would yield `TOLERANCE_MESH_LEGACY_f64`
# (grep the codebase for `TOLERANCE_.*_f64` if that happens).
REPLACEMENTS: list[tuple[str, str]] = [
    ("2.0e-6", "2.0 * TOLERANCE_MESH_LEGACY"),
    ("2.5e-7", "2.5 * TOLERANCE_ABS"),
    ("2.0e-7", "2.0 * TOLERANCE_ABS"),
    ("3.0e-7", "3.0 * TOLERANCE_ABS"),
    ("4.0e-7", "4.0 * TOLERANCE_ABS"),
    ("2.0e-5", "2.0 * TOLERANCE_RETRY_LADDER_MID"),
    ("1.0e-18", "TOLERANCE_FLOAT_ULTRA"),
    ("1.0e-15", "TOLERANCE_FLOAT_DEDUP"),
    ("1.0e-14", "TOLERANCE_FLOAT_LOOSE"),
    ("1.0e-12", "TOLERANCE_LEN_MIN"),
    ("1.0e-10", "TOLERANCE_LINEAR_ULTRA_STRICT"),
    ("1.0e-9", "TOLERANCE_COORD_SUB"),
    ("1.0e-8", "TOLERANCE_LINEAR_RELAX_8"),
    ("1.0e-7", "TOLERANCE_ABS"),
    ("1.0e-6", "TOLERANCE_MESH_LEGACY"),
    ("1.0e-5", "TOLERANCE_RETRY_LADDER_MID"),
    ("1.0e-4", "TOLERANCE_RETRY_LADDER_COARSE"),
    ("1.0e-3", "TOLERANCE_ADAPTIVE_MAX"),
    ("2e-6", "2.0 * TOLERANCE_MESH_LEGACY"),
    ("3e-7", "3.0 * TOLERANCE_ABS"),
    ("4e-7", "4.0 * TOLERANCE_ABS"),
    ("5e-4", "0.5 * TOLERANCE_RETRY_LADDER_COARSE"),
    ("5e-3", "50.0 * TOLERANCE_RETRY_LADDER_COARSE"),  # 5e-3 = 50 * 1e-4
    ("5e-6", "50.0 * TOLERANCE_ABS"),
    ("1e-20", "TOLERANCE_METRIC_SQ_NEAR_ZERO"),
    ("1e-30", "TOLERANCE_LEN_SQ_DIV_SAFE"),
    ("1e-24", "TOLERANCE_VEC_SQ_MIN"),
    ("1e-18", "TOLERANCE_FLOAT_ULTRA"),
    ("-1e-10", "-TOLERANCE_LINEAR_ULTRA_STRICT"),
    ("1e-15", "TOLERANCE_FLOAT_DEDUP"),
    ("1e-14", "TOLERANCE_FLOAT_LOOSE"),
    ("1e-12", "TOLERANCE_LEN_MIN"),
    ("1e-10", "TOLERANCE_LINEAR_ULTRA_STRICT"),
    ("1e-9", "TOLERANCE_COORD_SUB"),
    ("1e-8", "TOLERANCE_LINEAR_RELAX_8"),
    ("1e-7", "TOLERANCE_ABS"),
    ("2e-3", "2.0 * TOLERANCE_ADAPTIVE_MAX"),
    ("2e-2", "TOLERANCE_AXIS_CORNER_SLACK"),
    ("1e-3", "TOLERANCE_ADAPTIVE_MAX"),
    ("1e-4", "TOLERANCE_RETRY_LADDER_COARSE"),
    ("1e-5", "TOLERANCE_RETRY_LADDER_MID"),
    ("1e-6", "TOLERANCE_MESH_LEGACY"),
]


def process_file(path: Path) -> bool:
    if path.name == "tolerance.rs":
        return False
    raw = path.read_text(encoding="utf-8")
    out = raw
    for old, new in REPLACEMENTS:
        out = out.replace(old, new)
    if out == raw:
        return False
    path.write_text(out, encoding="utf-8")
    return True


def main() -> None:
    changed = 0
    for root in (SRC, TESTS):
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.rs")):
            if process_file(path):
                changed += 1
                print(path.relative_to(CRATE_ROOT.parent.parent.parent))
    print(f"updated {changed} files")


if __name__ == "__main__":
    main()
