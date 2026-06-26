with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find impl block start
impl_start = None
for i, line in enumerate(lines):
    if "impl<" in line and "PaveFiller" in line:
        impl_start = i
        break

if impl_start is not None:
    print(f"Impl block at line {impl_start+1}")
    depth = 0
    for i in range(impl_start, len(lines)):
        in_str = False
        for j, ch in enumerate(lines[i]):
            if ch == '"' and (j == 0 or lines[i][j-1] != '\\'):
                in_str = not in_str
            elif ch == '{' and not in_str:
                depth += 1
            elif ch == '}' and not in_str:
                depth -= 1
        if depth < 0:
            s = lines[i][:70]
            s = "".join(c if ord(c) < 128 else '?' for c in s)
            print(f"  Depth negative at line {i+1}: {s}")
            # Find where the imbalance is
            # Check the area around the opening
            print(f"  Impl started at line {impl_start+1}")
            break
    print(f"Final depth at end: {depth}")
