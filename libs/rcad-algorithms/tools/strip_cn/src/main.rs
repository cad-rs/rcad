use std::fs;
use std::path::Path;

fn main() {
    let path = Path::new("../../src/pave_filler/mod.rs");
    let content = fs::read_to_string(path).expect("read file");
    // Remove Chinese characters (CJK Unified Ideographs block U+4E00..U+9FFF)
    let cleaned: String = content
        .chars()
        .filter(|&c| !matches!(c as u32, 0x4E00..=0x9FFF | 0x3000..=0x303F | 0xFF00..=0xFFEF))
        .collect();
    fs::write(path, cleaned).expect("write file");
    println!("Chinese characters stripped from {}", path.display());
}
