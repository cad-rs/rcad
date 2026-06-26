with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Remove line 2554 (0-indexed: 2553) - the orphaned }
print(f"Line 2554: {lines[2553].rstrip()}")
del lines[2553]

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'w', encoding='utf-8') as f:
    f.writelines(lines)

print('Removed orphaned } at line 2554')

# Re-check brace balance
depth = 0
for i, line in enumerate(lines):
    for j, ch in enumerate(line):
        if ch == '"' and (j == 0 or line[j-1] != chr(92)):
            pass  # skip string content
    depth += line.count('{') - line.count('}')
    if depth < 0:
        s = ''.join(c if ord(c) < 128 else '?' for c in line[:70])
        print(f'Negative depth at line {i+1}: {s}')
        break

print(f'Final depth: {depth}')
