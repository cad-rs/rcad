//! Unit system — analogous to OCCT `Units` / `UnitsAPI`.
//!
//! Provides unit definitions, conversion factors, and a central conversion
//! API for length, angle, mass, time, and other physical dimensions.
//!
//! # Base units (SI)
//!
//! | Dimension | SI unit | Symbol |
//! |-----------|---------|--------|
//! | Length    | meter   | m      |
//! | Angle     | radian  | rad    |
//! | Mass      | kilogram| kg     |
//! | Time      | second  | s      |
//!
//! # Examples
//!
//! ```
//! use rcad_kernel::core::units::*;
//!
//! // Convert 1 inch to millimeters
//! let mm = to_si(1.0, "INCH").unwrap();        // → 0.0254 (meters)
//! let inch = from_si(mm, "INCH").unwrap();     // → 1.0
//!
//! // Convert between arbitrary units
//! let mm_per_inch = convert(1.0, "INCH", "MM").unwrap(); // → 25.4
//!
//! // Set current unit system
//! set_current_unit_system(UnitSystem::SI);
//! assert_eq!(current_length_unit(), "MM");
//! ```

use std::collections::HashMap;
use std::sync::Mutex;

// ============================================================================
// Unit dimensions
// ============================================================================

/// Physical dimension of a unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum UnitDimension {
    Length,
    Angle,
    Mass,
    Time,
    Force,
    Pressure,
    Area,
    Volume,
    Velocity,
    Temperature,
}

// ============================================================================
// Unit system
// ============================================================================

/// Predefined unit systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnitSystem {
    /// SI: meter, radian, kilogram, second.
    SI,
    /// Millimeter-based: millimeter, radian, kilogram, second.
    /// Common in mechanical CAD.
    MM,
    /// Inch-based: inch, radian, pound, second.
    /// Common in imperial CAD workflows.
    Inch,
}

impl UnitSystem {
    /// The length unit name for this system.
    pub fn length_unit(self) -> &'static str {
        match self {
            UnitSystem::SI => "M",
            UnitSystem::MM => "MM",
            UnitSystem::Inch => "INCH",
        }
    }

    /// The angle unit name for this system.
    pub fn angle_unit(self) -> &'static str {
        "RAD"
    }

    /// The mass unit name for this system.
    pub fn mass_unit(self) -> &'static str {
        match self {
            UnitSystem::SI | UnitSystem::MM => "KG",
            UnitSystem::Inch => "LB",
        }
    }
}

// ============================================================================
// Conversion factors to SI
// ============================================================================
//
// Each entry: value_in_SI = value * factor
// For example, to_si_factor("MM") = 0.001 (1 mm = 0.001 m)

fn build_conversion_table() -> HashMap<&'static str, (UnitDimension, f64)> {
    let mut t: HashMap<&'static str, (UnitDimension, f64)> = HashMap::new();

    // ── Length ──────────────────────────────────────────────────────────
    t.insert("M", (UnitDimension::Length, 1.0));
    t.insert("METRE", (UnitDimension::Length, 1.0));
    t.insert("METER", (UnitDimension::Length, 1.0));
    t.insert("MM", (UnitDimension::Length, 1e-3));
    t.insert("MILLIMETRE", (UnitDimension::Length, 1e-3));
    t.insert("MILLIMETER", (UnitDimension::Length, 1e-3));
    t.insert("CM", (UnitDimension::Length, 1e-2));
    t.insert("CENTIMETRE", (UnitDimension::Length, 1e-2));
    t.insert("KM", (UnitDimension::Length, 1e3));
    t.insert("INCH", (UnitDimension::Length, 0.0254));
    t.insert("FT", (UnitDimension::Length, 0.3048));
    t.insert("FOOT", (UnitDimension::Length, 0.3048));
    t.insert("MI", (UnitDimension::Length, 1609.344));
    t.insert("MICRON", (UnitDimension::Length, 1e-6));

    // ── Angle ───────────────────────────────────────────────────────────
    t.insert("RAD", (UnitDimension::Angle, 1.0));
    t.insert("RADIAN", (UnitDimension::Angle, 1.0));
    t.insert("DEG", (UnitDimension::Angle, std::f64::consts::PI / 180.0));
    t.insert("DEGREE", (UnitDimension::Angle, std::f64::consts::PI / 180.0));

    // ── Mass ────────────────────────────────────────────────────────────
    t.insert("KG", (UnitDimension::Mass, 1.0));
    t.insert("KILOGRAM", (UnitDimension::Mass, 1.0));
    t.insert("G", (UnitDimension::Mass, 1e-3));
    t.insert("GRAM", (UnitDimension::Mass, 1e-3));
    t.insert("LB", (UnitDimension::Mass, 0.45359237));
    t.insert("POUND", (UnitDimension::Mass, 0.45359237));
    t.insert("TONNE", (UnitDimension::Mass, 1000.0));

    // ── Time ────────────────────────────────────────────────────────────
    t.insert("S", (UnitDimension::Time, 1.0));
    t.insert("SECOND", (UnitDimension::Time, 1.0));
    t.insert("MIN", (UnitDimension::Time, 60.0));
    t.insert("MINUTE", (UnitDimension::Time, 60.0));
    t.insert("HR", (UnitDimension::Time, 3600.0));
    t.insert("HOUR", (UnitDimension::Time, 3600.0));

    // ── Force ───────────────────────────────────────────────────────────
    t.insert("N", (UnitDimension::Force, 1.0));
    t.insert("NEWTON", (UnitDimension::Force, 1.0));
    t.insert("LBF", (UnitDimension::Force, 4.4482216152605));

    t
}

fn conversion_table() -> &'static HashMap<&'static str, (UnitDimension, f64)> {
    static TABLE: std::sync::LazyLock<HashMap<&'static str, (UnitDimension, f64)>> =
        std::sync::LazyLock::new(build_conversion_table);
    &TABLE
}

// ============================================================================
// Core conversion functions
// ============================================================================

/// Look up the conversion factor to SI for a named unit.
///
/// Returns `(dimension, factor)` where `value_in_SI = value * factor`.
pub fn to_si_factor(unit: &str) -> Option<(UnitDimension, f64)> {
    let upper = unit.to_uppercase();
    conversion_table().get(upper.as_str()).copied()
}

/// Convert a value from a named unit to SI.
///
/// # Example
/// ```
/// let meters = to_si(25.4, "MM").unwrap();
/// assert!((meters - 0.0254).abs() < 1e-12);
/// ```
pub fn to_si(value: f64, unit: &str) -> Option<f64> {
    to_si_factor(unit).map(|(_, factor)| value * factor)
}

/// Convert a value from SI to a named unit.
///
/// # Example
/// ```
/// let mm = from_si(0.0254, "MM").unwrap();
/// assert!((mm - 25.4).abs() < 1e-12);
/// ```
pub fn from_si(value: f64, unit: &str) -> Option<f64> {
    to_si_factor(unit).map(|(_, factor)| value / factor)
}

/// Convert a value between two arbitrary units.
///
/// Returns `None` if either unit name is unknown or the dimensions differ.
///
/// # Example
/// ```
/// let mm_per_inch = convert(1.0, "INCH", "MM").unwrap();
/// assert!((mm_per_inch - 25.4).abs() < 1e-12);
/// ```
pub fn convert(value: f64, from: &str, to: &str) -> Option<f64> {
    let table = conversion_table();
    let upper_from = from.to_uppercase();
    let upper_to = to.to_uppercase();

    let (dim_from, factor_from) = table.get(upper_from.as_str())?;
    let (dim_to, factor_to) = table.get(upper_to.as_str())?;

    if dim_from != dim_to {
        return None;
    }

    // value (from_units) → value * factor_from (SI) → (value * factor_from) / factor_to (to_units)
    Some(value * factor_from / factor_to)
}

/// Check whether two units represent the same dimension.
pub fn same_dimension(a: &str, b: &str) -> bool {
    let table = conversion_table();
    let upper_a = a.to_uppercase();
    let upper_b = b.to_uppercase();

    match (table.get(upper_a.as_str()), table.get(upper_b.as_str())) {
        (Some((dim_a, _)), Some((dim_b, _))) => dim_a == dim_b,
        _ => false,
    }
}

// ============================================================================
// Current unit system (UnitsAPI style)
// ============================================================================

static CURRENT_SYSTEM: Mutex<UnitSystem> = Mutex::new(UnitSystem::MM);

/// Set the current unit system.
pub fn set_current_unit_system(system: UnitSystem) {
    *CURRENT_SYSTEM.lock().unwrap() = system;
}

/// Get the current unit system.
pub fn current_unit_system() -> UnitSystem {
    *CURRENT_SYSTEM.lock().unwrap()
}

/// Get the current length unit name.
pub fn current_length_unit() -> &'static str {
    current_unit_system().length_unit()
}

/// Convert a value from the current system's length unit to SI (meters).
pub fn to_si_length(value: f64) -> f64 {
    let unit = current_length_unit();
    to_si(value, unit).expect("current length unit should be in conversion table")
}

/// Convert a value from SI (meters) to the current system's length unit.
pub fn from_si_length(value: f64) -> f64 {
    let unit = current_length_unit();
    from_si(value, unit).expect("current length unit should be in conversion table")
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mm_to_m() {
        let m = to_si(1000.0, "MM").unwrap();
        assert!((m - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_inch_to_mm() {
        let mm = convert(1.0, "INCH", "MM").unwrap();
        assert!((mm - 25.4).abs() < 1e-12);
    }

    #[test]
    fn test_mm_to_inch() {
        let inch = convert(25.4, "MM", "INCH").unwrap();
        assert!((inch - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_degree_to_radian() {
        let rad = to_si(180.0, "DEG").unwrap();
        assert!((rad - std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn test_roundtrip() {
        let values = [1.0, 25.4, 1000.0, 0.001];
        let units = ["MM", "INCH", "M", "CM", "FT"];
        for &v in &values {
            for &u in &units {
                let si = to_si(v, u).unwrap();
                let back = from_si(si, u).unwrap();
                assert!(
                    (back - v).abs() < 1e-12,
                    "roundtrip failed for {} in {}: {} → {} → {}",
                    v,
                    u,
                    v,
                    si,
                    back
                );
            }
        }
    }

    #[test]
    fn test_dimension_mismatch_returns_none() {
        assert!(convert(1.0, "MM", "DEG").is_none());
        assert!(convert(1.0, "KG", "S").is_none());
    }

    #[test]
    fn test_unknown_unit_returns_none() {
        assert!(to_si(1.0, "FURLONG").is_none());
        assert!(from_si(1.0, "FURLONG").is_none());
    }

    #[test]
    fn test_case_insensitive() {
        let a = to_si(1.0, "mm").unwrap();
        let b = to_si(1.0, "MM").unwrap();
        assert!((a - b).abs() < 1e-12);
    }

    #[test]
    fn test_current_system_default() {
        assert_eq!(current_unit_system(), UnitSystem::MM);
        set_current_unit_system(UnitSystem::SI);
        assert_eq!(current_unit_system(), UnitSystem::SI);
        set_current_unit_system(UnitSystem::MM); // restore
    }

    #[test]
    fn test_to_si_length() {
        set_current_unit_system(UnitSystem::MM);
        let m = to_si_length(1000.0);
        assert!((m - 1.0).abs() < 1e-12);

        set_current_unit_system(UnitSystem::SI);
        let m = to_si_length(1.0);
        assert!((m - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_same_dimension() {
        assert!(same_dimension("MM", "INCH"));
        assert!(same_dimension("M", "FT"));
        assert!(!same_dimension("MM", "DEG"));
        assert!(!same_dimension("MM", "unknown"));
    }

    #[test]
    fn test_all_known_units_roundtrip() {
        let known = ["M", "MM", "CM", "INCH", "FT", "RAD", "DEG", "KG", "G", "LB", "S", "MIN"];
        for &unit in &known {
            let si = to_si(1.0, unit);
            assert!(si.is_some(), "unit {} not found", unit);
            let back = from_si(si.unwrap(), unit);
            assert!((back.unwrap() - 1.0).abs() < 1e-12, "unit {}", unit);
        }
    }
}
