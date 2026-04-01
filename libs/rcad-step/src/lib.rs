use rcad_kernel::BRep;
use std::path::Path;

pub struct StepReader;

impl StepReader {
    pub fn read_file<P: AsRef<Path>>(_path: P) -> Result<BRep, String> {
        // Placeholder for actual STEP parsing
        // For now, return a mock box to unblock UI development
        Ok(BRep::create_box(10.0, 10.0, 10.0))
    }

    pub fn parse_string(content: &str) -> Result<BRep, String> {
        if content.contains("ISO-10303-21") {
            Ok(BRep::create_box(1.0, 1.0, 1.0))
        } else {
            Err("Invalid STEP file format".to_string())
        }
    }
}
