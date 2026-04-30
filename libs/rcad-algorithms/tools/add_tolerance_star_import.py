"""Insert `use crate::tolerance::*;` in modules that reference TOLERANCE_* but lack a tolerance import."""
from __future__ import annotations

import re
from pathlib import Path

CRATE_ROOT = Path(__file__).resolve().parents[1]
SRC = CRATE_ROOT / "src"
TESTS = CRATE_ROOT / "tests"


def needs_import(text: str, is_integration: bool) -> bool:
    if "TOLERANCE_" not in text:
        return False
    if re.search(r"^\s*use\s+crate::tolerance::", text, re.MULTILINE):
        return False
    if is_integration and re.search(
        r"^\s*use\s+rcad_algorithms::tolerance::", text, re.MULTILINE
    ):
        return False
    return True


def insert_import(text: str, is_integration: bool) -> str:
    lines = text.splitlines(keepends=True)
    idx = 0
    n = len(lines)
    while idx < n and (
        lines[idx].startswith("#!")
        or lines[idx].startswith("//!")
        or lines[idx].startswith("///")
        or lines[idx].strip() == ""
    ):
        idx += 1
    import_line = (
        "use rcad_algorithms::tolerance::*;\n"
        if is_integration
        else "use crate::tolerance::*;\n"
    )
    lines.insert(idx, import_line)
    return "".join(lines)


def process(path: Path) -> bool:
    if path.name == "tolerance.rs":
        return False
    is_integration = path.parent == TESTS
    raw = path.read_text(encoding="utf-8")
    if not needs_import(raw, is_integration):
        return False
    path.write_text(insert_import(raw, is_integration), encoding="utf-8")
    return True


def main() -> None:
    n = 0
    for root in (SRC, TESTS):
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.rs")):
            if process(path):
                n += 1
                print(path.relative_to(CRATE_ROOT))
    print(f"inserted import in {n} files")


if __name__ == "__main__":
    main()
