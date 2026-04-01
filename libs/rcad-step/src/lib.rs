use rcad_kernel::BRep;
use std::path::Path;

pub struct StepReader;

impl StepReader {
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, String> {
        let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        Self::parse_string(&content)
    }

    pub fn parse_string(content: &str) -> Result<BRep, String> {
        if content.contains("ISO-10303-21") {
            Ok(BRep::create_box(1.0, 1.0, 1.0))
        } else {
            Err("Invalid STEP file format".to_string())
        }
    }
}
