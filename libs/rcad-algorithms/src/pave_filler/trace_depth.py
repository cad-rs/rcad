with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'r', encoding='utf-8') as f:
    lines = f.readlines()

impl_start = None
for i, line in enumerate(lines):
    if "impl<" in line and "PaveFiller" in line:
        impl_start = i
        break

if impl_start is None:
    print("ERROR: impl block not found")
    exit(1)

print(f"Impl block at line {impl_start+1}")
depth = 0
for i in range(impl_start, len(lines)):
    depth += lines[i].count("{") - lines[i].count("}")
    if i+1 == 2367:
        print(f"  At line {i+1} (error location): depth = {depth}")
    if depth < 0:
        print(f"  Depth negative at line {i+1}: {lines[i][:50]}")
        # Find where this } came from
        break

print(f"Final depth at end of file: {depth}")
