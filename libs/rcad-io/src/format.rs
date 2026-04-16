/// Supported layout file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LayoutFormat {
    /// GDSII format (.gds, .gds2)
    Gds,
    /// OASIS format (.oas, .oasis)
    Oasis,
}

impl LayoutFormat {
    /// Returns the default file extension for this format.
    pub fn extension(&self) -> &'static str {
        match self {
            LayoutFormat::Gds => "gds",
            LayoutFormat::Oasis => "oas",
        }
    }

    /// Returns all valid extensions for this format.
    pub fn extensions(&self) -> &'static [&'static str] {
        match self {
            LayoutFormat::Gds => &["gds", "gds2"],
            LayoutFormat::Oasis => &["oas", "oasis"],
        }
    }

    /// Detects format from file extension.
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_lowercase().as_str() {
            "gds" | "gds2" => Some(LayoutFormat::Gds),
            "oas" | "oasis" => Some(LayoutFormat::Oasis),
            _ => None,
        }
    }

    /// Returns the MIME type for this format.
    pub fn mime_type(&self) -> &'static str {
        match self {
            LayoutFormat::Gds => "application/x-gds2",
            LayoutFormat::Oasis => "application/x-oasis",
        }
    }

    /// Returns a human-readable name for this format.
    pub fn name(&self) -> &'static str {
        match self {
            LayoutFormat::Gds => "GDSII",
            LayoutFormat::Oasis => "OASIS",
        }
    }
}

impl std::fmt::Display for LayoutFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}
