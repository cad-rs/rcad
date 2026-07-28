#!/usr/bin/env python3
"""Convert ShapeRef -> Shape in topods.rs with proper UTF-8 handling."""

import re

FILE = r"C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-kernel\src\topods.rs"

# Read as bytes, then decode as UTF-8 (strict)
with open(FILE, "rb") as f:
    raw = f.read()
content = raw.decode("utf-8")

original = content

# Step 1: Remove the SYNTH_PTR_ID line
content = re.sub(r'/\* Sentinel ptr_id value.*?\n', '', content, flags=re.DOTALL)

# Step 2: Remove pub struct ShapeRef and all impl blocks for ShapeRef
# Remove the struct definition
content = re.sub(
    r'/\*\* TopoDS_Shape equivalent.*?pub struct ShapeRef \{.*?\n\}',
    '', content, flags=re.DOTALL
)

# Remove impl PartialEq/Eq/Hash/Ord/Serialize/Deserialize for ShapeRef
content = re.sub(
    r'impl (PartialEq|Eq|std::hash::Hash|PartialOrd|Ord|Serialize) for ShapeRef \{.*?\n\}',
    '', content, flags=re.DOTALL
)

content = re.sub(
    r'impl<.*?> Deserialize<.*?> for ShapeRef \{.*?\n\}',
    '', content, flags=re.DOTALL
)

# Remove impl ShapeRef block
content = re.sub(
    r'impl ShapeRef \{.*?\n\}',
    '', content, flags=re.DOTALL
)

# Step 3: Replace ShapeRef method calls with Shape
content = content.replace('ShapeRef::NULL', 'Shape::null()')
content = content.replace('ShapeRef::synthetic_with_orientation(', 'Shape::synthetic(')
content = content.replace('ShapeRef::synthetic_with_location(', 'Shape::synthetic_with_location(')
content = content.replace('ShapeRef::synthetic(', 'Shape::synthetic(')

# Step 4: Replace type references
content = re.sub(r'->\s*ShapeRef\b', '-> Shape', content)
content = re.sub(r'(?<!::)\bShapeRef\b(?!::)', 'Shape', content)

# Step 5: Replace string literals
content = content.replace('"ShapeRef', '"Shape')

# Step 6: Fix struct literals - ShapeRef { ptr_id, ... } -> Shape { data: ..., ... }
# The old struct had ptr_id field, new one has data field.
# Replace "Shape {" that has ptr_id: in it
# Pattern: Shape { ptr_id: ... -> need to fix
# These need manual handling since they vary per context.

if content != original:
    with open(FILE, "wb") as f:
        f.write(content.encode("utf-8"))
    print("File updated successfully")
else:
    print("No changes made")

# Count remaining ShapeRef
remaining = content.count('ShapeRef')
print(f"Remaining 'ShapeRef': {remaining}")
