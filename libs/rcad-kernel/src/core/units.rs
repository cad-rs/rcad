//! Unit system — analogous to OCCT `Units` / `UnitsAPI`.
//!
//! OCCT defines three co-existing unit systems:
//!
//! - **SI System**: standard international (meter, kilogram, second, radian, …)
//! - **Local System (LS)**: either SI or MDTV.
//!   MDTV = SI with millimeter as the length base (mm, kg, s). This is the
//!   historical Open CASCADE default.
//! - **Current System**: user-customizable per quantity (e.g. length in inches,
//!   angle in degrees) via `set_current_unit()`.
//!
//! The unit dictionary is built-in (hard-coded, not loaded from a file).
//!
//! # Architecture
//!
//! | OCCT class       | rcad equivalent              | Role                          |
//! |------------------|------------------------------|-------------------------------|
//! | `Units_Dimensions` | `Dimensions`                | Physical dimension (M, L, T…) |
//! | `Units_Unit`      | `Unit`                      | Single unit with SI factor    |
//! | `Units_Quantity`  | `Quantity`                  | A physical quantity kind      |
//! | `Units_UnitsDictionary` | `UnitsDictionary`    | Registry of quantities+units  |
//! | `Units`           | `Units` (static)            | `to_si()`, `from_si()`        |
//! | `UnitsAPI`        | `UnitsAPI` (static)         | Three-system conversion API   |

use std::collections::HashMap;
use std::sync::Mutex;

// ============================================================================
// Units_Dimensions
// ============================================================================

/// Physical dimensions expressed as exponents of 9 base quantities.
///
/// OCCT: `Units_Dimensions` (M, L, T, I, Θ, N, J, plane angle, solid angle).
#[derive(Debug, Clone, PartialEq)]
pub struct Dimensions {
    pub mass: f64,
    pub length: f64,
    pub time: f64,
    pub electric_current: f64,
    pub thermodynamic_temperature: f64,
    pub amount_of_substance: f64,
    pub luminous_intensity: f64,
    pub plane_angle: f64,
    pub solid_angle: f64,
}

impl Dimensions {
    /// OCCT: `Units_Dimensions(M, L, T, I, Θ, N, J, plane_angle, solid_angle)`.
    pub fn new(
        mass: f64,
        length: f64,
        time: f64,
        electric_current: f64,
        thermodynamic_temperature: f64,
        amount_of_substance: f64,
        luminous_intensity: f64,
        plane_angle: f64,
        solid_angle: f64,
    ) -> Self {
        Self {
            mass,
            length,
            time,
            electric_current,
            thermodynamic_temperature,
            amount_of_substance,
            luminous_intensity,
            plane_angle,
            solid_angle,
        }
    }

    /// OCCT: `ALength()` — dimension with Length = 1.
    pub fn length() -> Self {
        Self::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }
    /// OCCT: `AMass()`.
    pub fn mass() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }
    /// OCCT: `ATime()`.
    pub fn time() -> Self {
        Self::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }
    /// OCCT: `APlaneAngle()`.
    pub fn plane_angle() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    }
    /// Dimensionless.
    pub fn less() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }
    /// OCCT: `Multiply(adim)`.
    pub fn multiply(&self, other: &Dimensions) -> Self {
        Self {
            mass: self.mass + other.mass,
            length: self.length + other.length,
            time: self.time + other.time,
            electric_current: self.electric_current + other.electric_current,
            thermodynamic_temperature: self.thermodynamic_temperature + other.thermodynamic_temperature,
            amount_of_substance: self.amount_of_substance + other.amount_of_substance,
            luminous_intensity: self.luminous_intensity + other.luminous_intensity,
            plane_angle: self.plane_angle + other.plane_angle,
            solid_angle: self.solid_angle + other.solid_angle,
        }
    }
    /// OCCT: `Divide(adim)`.
    pub fn divide(&self, other: &Dimensions) -> Self {
        Self {
            mass: self.mass - other.mass,
            length: self.length - other.length,
            time: self.time - other.time,
            electric_current: self.electric_current - other.electric_current,
            thermodynamic_temperature: self.thermodynamic_temperature - other.thermodynamic_temperature,
            amount_of_substance: self.amount_of_substance - other.amount_of_substance,
            luminous_intensity: self.luminous_intensity - other.luminous_intensity,
            plane_angle: self.plane_angle - other.plane_angle,
            solid_angle: self.solid_angle - other.solid_angle,
        }
    }
    /// OCCT: `Power(exp)`.
    pub fn power(&self, exp: f64) -> Self {
        Self {
            mass: self.mass * exp,
            length: self.length * exp,
            time: self.time * exp,
            electric_current: self.electric_current * exp,
            thermodynamic_temperature: self.thermodynamic_temperature * exp,
            amount_of_substance: self.amount_of_substance * exp,
            luminous_intensity: self.luminous_intensity * exp,
            plane_angle: self.plane_angle * exp,
            solid_angle: self.solid_angle * exp,
        }
    }
    /// OCCT: `IsEqual(adim)`.
    pub fn is_equal(&self, other: &Dimensions) -> bool {
        const TOL: f64 = 1e-12;
        (self.mass - other.mass).abs() < TOL
            && (self.length - other.length).abs() < TOL
            && (self.time - other.time).abs() < TOL
            && (self.electric_current - other.electric_current).abs() < TOL
            && (self.thermodynamic_temperature - other.thermodynamic_temperature).abs() < TOL
            && (self.amount_of_substance - other.amount_of_substance).abs() < TOL
            && (self.luminous_intensity - other.luminous_intensity).abs() < TOL
            && (self.plane_angle - other.plane_angle).abs() < TOL
            && (self.solid_angle - other.solid_angle).abs() < TOL
    }
}

// ============================================================================
// Units_Unit
// ============================================================================

/// A single unit with its name, symbol(s), and SI conversion factor.
///
/// OCCT: `Units_Unit(name, symbol, value, quantity)`.
#[derive(Debug, Clone)]
pub struct Unit {
    pub name: String,
    pub symbols: Vec<String>,
    /// Conversion factor to SI: `value_in_SI = value * factor`.
    pub factor: f64,
    pub quantity_name: String,
}

impl Unit {
    pub fn new(name: &str, symbol: &str, factor: f64, quantity: &str) -> Self {
        Self {
            name: name.to_string(),
            symbols: vec![symbol.to_string()],
            factor,
            quantity_name: quantity.to_string(),
        }
    }
}

// ============================================================================
// Units_Quantity
// ============================================================================

/// A physical quantity (e.g. "Length", "Mass") grouping its possible units.
///
/// OCCT: `Units_Quantity(name, dimensions, units_sequence)`.
#[derive(Debug, Clone)]
pub struct Quantity {
    pub name: String,
    pub dimensions: Dimensions,
    pub units: Vec<Unit>,
}

impl Quantity {
    pub fn new(name: &str, dimensions: Dimensions, units: Vec<Unit>) -> Self {
        Self {
            name: name.to_string(),
            dimensions,
            units,
        }
    }
}

// ============================================================================
// Units_UnitsDictionary
// ============================================================================

/// Registry of all known physical quantities and their units.
///
/// OCCT: `Units_UnitsDictionary`.
/// Built-in (no file loading).
#[derive(Debug, Clone)]
pub struct UnitsDictionary {
    quantities: Vec<Quantity>,
}

impl UnitsDictionary {
    /// OCCT: `Creates()` — build the dictionary.
    pub fn new() -> Self {
        Self::create()
    }

    /// Build the built-in unit dictionary.
    fn create() -> Self {
        let mut dict = Self {
            quantities: Vec::new(),
        };

        // ── Length ────────────────────────────────────────────────────
        dict.quantities.push(Quantity::new(
            "Length",
            Dimensions::length(),
            vec![
                Unit::new("meter", "m", 1.0, "Length"),
                Unit::new("millimeter", "mm", 1e-3, "Length"),
                Unit::new("centimeter", "cm", 1e-2, "Length"),
                Unit::new("kilometer", "km", 1e3, "Length"),
                Unit::new("inch", "in", 0.0254, "Length"),
                Unit::new("foot", "ft", 0.3048, "Length"),
                Unit::new("mile", "mi", 1609.344, "Length"),
                Unit::new("micron", "µm", 1e-6, "Length"),
            ],
        ));

        // ── Mass ──────────────────────────────────────────────────────
        dict.quantities.push(Quantity::new(
            "Mass",
            Dimensions::mass(),
            vec![
                Unit::new("kilogram", "kg", 1.0, "Mass"),
                Unit::new("gram", "g", 1e-3, "Mass"),
                Unit::new("pound", "lb", 0.45359237, "Mass"),
                Unit::new("tonne", "t", 1000.0, "Mass"),
            ],
        ));

        // ── Time ──────────────────────────────────────────────────────
        dict.quantities.push(Quantity::new(
            "Time",
            Dimensions::time(),
            vec![
                Unit::new("second", "s", 1.0, "Time"),
                Unit::new("minute", "min", 60.0, "Time"),
                Unit::new("hour", "h", 3600.0, "Time"),
            ],
        ));

        // ── PlaneAngle ────────────────────────────────────────────────
        dict.quantities.push(Quantity::new(
            "PlaneAngle",
            Dimensions::plane_angle(),
            vec![
                Unit::new("radian", "rad", 1.0, "PlaneAngle"),
                Unit::new("degree", "deg", std::f64::consts::PI / 180.0, "PlaneAngle"),
            ],
        ));

        // ── SolidAngle ────────────────────────────────────────────────
        dict.quantities.push(Quantity::new(
            "SolidAngle",
            Dimensions::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
            vec![
                Unit::new("steradian", "sr", 1.0, "SolidAngle"),
            ],
        ));

        // ── Velocity ──────────────────────────────────────────────────
        let vel_dim = Dimensions::length().multiply(&Dimensions::time().power(-1.0));
        dict.quantities.push(Quantity::new(
            "Velocity",
            vel_dim,
            vec![
                Unit::new("meter per second", "m/s", 1.0, "Velocity"),
                Unit::new("kilometer per hour", "km/h", 1.0 / 3.6, "Velocity"),
            ],
        ));

        // ── Force ─────────────────────────────────────────────────────
        let force_dim = Dimensions::mass()
            .multiply(&Dimensions::length())
            .multiply(&Dimensions::time().power(-2.0));
        let press_dim = force_dim.clone().multiply(&Dimensions::length().power(-2.0));
        dict.quantities.push(Quantity::new(
            "Force",
            force_dim,
            vec![
                Unit::new("newton", "N", 1.0, "Force"),
                Unit::new("pound-force", "lbf", 4.4482216152605, "Force"),
            ],
        ));

        // ── Pressure ──────────────────────────────────────────────────
        // press_dim computed above from the cloned force_dim
        dict.quantities.push(Quantity::new(
            "Pressure",
            press_dim,
            vec![
                Unit::new("pascal", "Pa", 1.0, "Pressure"),
                Unit::new("bar", "bar", 1e5, "Pressure"),
            ],
        ));

        // ── Area ──────────────────────────────────────────────────────
        let area_dim = Dimensions::length().power(2.0);
        dict.quantities.push(Quantity::new(
            "Area",
            area_dim,
            vec![
                Unit::new("square meter", "m2", 1.0, "Area"),
                Unit::new("square millimeter", "mm2", 1e-6, "Area"),
                Unit::new("square inch", "in2", 0.0254 * 0.0254, "Area"),
            ],
        ));

        // ── Volume ────────────────────────────────────────────────────
        let vol_dim = Dimensions::length().power(3.0);
        dict.quantities.push(Quantity::new(
            "Volume",
            vol_dim,
            vec![
                Unit::new("cubic meter", "m3", 1.0, "Volume"),
                Unit::new("cubic millimeter", "mm3", 1e-9, "Volume"),
                Unit::new("liter", "L", 1e-3, "Volume"),
            ],
        ));

        dict
    }

    /// OCCT: `Sequence()` — all quantities.
    pub fn quantities(&self) -> &[Quantity] {
        &self.quantities
    }

    /// Find a quantity by name (case-insensitive).
    pub fn find_quantity(&self, name: &str) -> Option<&Quantity> {
        let upper = name.to_uppercase();
        self.quantities.iter().find(|q| q.name.to_uppercase() == upper)
    }

    /// Find a unit by symbol (case-insensitive, across all quantities).
    pub fn find_unit(&self, symbol: &str) -> Option<&Unit> {
        let upper = symbol.to_uppercase();
        for q in &self.quantities {
            for u in &q.units {
                if u.symbols.iter().any(|s| s.to_uppercase() == upper) {
                    return Some(u);
                }
                if u.name.to_uppercase() == upper {
                    return Some(u);
                }
            }
        }
        None
    }

    /// OCCT: `ActiveUnit(quantity)` — the preferred unit for display.
    pub fn active_unit(&self, quantity: &str) -> Option<&str> {
        self.find_quantity(quantity)
            .and_then(|q| q.units.first())
            .map(|u| u.symbols[0].as_str())
    }
}

// ============================================================================
// Global dictionary singleton
// ============================================================================

/// OCCT: `Units::DictionaryOfUnits()`.
pub fn dictionary_of_units() -> &'static UnitsDictionary {
    static DICT: std::sync::LazyLock<UnitsDictionary> =
        std::sync::LazyLock::new(UnitsDictionary::new);
    &DICT
}

// ============================================================================
// Units::ToSI / FromSI / Convert
// ============================================================================

/// OCCT: `Units::ToSI(value, unit_name)`.
/// Convert from the given unit to SI.
pub fn to_si(value: f64, unit: &str) -> Option<f64> {
    dictionary_of_units()
        .find_unit(unit)
        .map(|u| value * u.factor)
}

/// OCCT: `Units::FromSI(value, unit_name)`.
/// Convert from SI to the given unit.
pub fn from_si(value: f64, unit: &str) -> Option<f64> {
    dictionary_of_units()
        .find_unit(unit)
        .map(|u| value / u.factor)
}

/// OCCT: `Units::Convert(value, from_unit, to_unit)`.
/// Convert between two arbitrary units.
///
/// Returns `None` when either unit is unknown or the two units belong to
/// different quantities (dimension mismatch) — OCCT `Units_Measurement::Convert`
/// refuses such a conversion.
pub fn convert(value: f64, from_unit: &str, to_unit: &str) -> Option<f64> {
    let dict = dictionary_of_units();
    let from = dict.find_unit(from_unit)?;
    let to = dict.find_unit(to_unit)?;
    if from.quantity_name != to.quantity_name {
        return None;
    }
    Some(value * from.factor / to.factor)
}

// ============================================================================
// Local System (LS): SI or MDTV
// ============================================================================

/// OCCT: `UnitsAPI_SystemUnits` — SI, MDTV (mm-based), or SI-by-default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SystemUnits {
    /// Standard international (meter, kilogram, second, …).
    SI,
    /// SI with millimeter for length and all derived units.
    /// This is the CAD-default system.
    MDTV,
}

static LOCAL_SYSTEM: Mutex<SystemUnits> = Mutex::new(SystemUnits::MDTV);

/// OCCT: `UnitsAPI::SetLocalSystem(system)`.
pub fn set_local_system(system: SystemUnits) {
    *LOCAL_SYSTEM.lock().unwrap() = system;
}

/// OCCT: `UnitsAPI::LocalSystem()`.
pub fn local_system() -> SystemUnits {
    *LOCAL_SYSTEM.lock().unwrap()
}

// ============================================================================
// Current System: per-quantity unit customization
// ============================================================================

fn current_units() -> &'static Mutex<HashMap<String, String>> {
    static MAP: std::sync::LazyLock<Mutex<HashMap<String, String>>> =
        std::sync::LazyLock::new(|| Mutex::new(HashMap::new()));
    &MAP
}

/// OCCT: `UnitsAPI::SetCurrentUnit(quantity_name, unit_name)`.
pub fn set_current_unit(quantity: &str, unit: &str) {
    let mut map = current_units().lock().unwrap();
    map.insert(quantity.to_string(), unit.to_string());
}

/// OCCT: `UnitsAPI::CurrentUnit(quantity_name)`.
/// Returns the current unit for the given quantity, or the default LS unit.
pub fn current_unit(quantity: &str) -> String {
    let map = current_units().lock().unwrap();
    if let Some(u) = map.get(quantity) {
        return u.clone();
    }
    // Fall back to LS default
    default_unit_for_quantity(quantity, local_system())
}

fn default_unit_for_quantity(quantity: &str, ls: SystemUnits) -> String {
    let dict = dictionary_of_units();
    if let Some(q) = dict.find_quantity(quantity) {
        if let Some(first) = q.units.first() {
            // For MDTV, prefer millimeter for Length-related quantities
            if ls == SystemUnits::MDTV && quantity.eq_ignore_ascii_case("Length") {
                return "mm".to_string();
            }
            return first.symbols[0].clone();
        }
    }
    // Fallback unknown quantities to SI unit
    quantity.to_string()
}

// ============================================================================
// UnitsAPI: three-system conversion
// ============================================================================

/// Converts the local system units value to SI.
/// OCCT: `UnitsAPI::LSToSI(data, quantity)`.
pub fn ls_to_si(data: f64, quantity: &str) -> Option<f64> {
    let unit = default_unit_for_quantity(quantity, local_system());
    to_si(data, &unit)
}

/// Converts SI to local system units.
/// OCCT: `UnitsAPI::SIToLS(data, quantity)`.
pub fn si_to_ls(data: f64, quantity: &str) -> Option<f64> {
    let unit = default_unit_for_quantity(quantity, local_system());
    from_si(data, &unit)
}

/// Converts from Current System to SI.
/// OCCT: `UnitsAPI::CurrentToSI(data, quantity)`.
pub fn current_to_si(data: f64, quantity: &str) -> Option<f64> {
    let unit = current_unit(quantity);
    to_si(data, &unit)
}

/// Converts from SI to Current System.
/// OCCT: `UnitsAPI::CurrentFromSI(data, quantity)`.
pub fn current_from_si(data: f64, quantity: &str) -> Option<f64> {
    let unit = current_unit(quantity);
    from_si(data, &unit)
}

/// Converts from any unit to LS.
/// OCCT: `UnitsAPI::AnyToLS(data, unit_name)`.
pub fn any_to_ls(data: f64, unit_name: &str) -> Option<f64> {
    let ls_unit = default_unit_for_quantity(
        dictionary_of_units().find_unit(unit_name)?.quantity_name.as_str(),
        local_system(),
    );
    let from = to_si(data, unit_name)?;
    from_si(from, &ls_unit)
}

/// Converts from any unit to SI.
/// OCCT: `UnitsAPI::AnyToSI(data, unit_name)`.
pub fn any_to_si(data: f64, unit_name: &str) -> Option<f64> {
    to_si(data, unit_name)
}

/// Converts from SI to any unit.
/// OCCT: `UnitsAPI::AnyFromSI(data, unit_name)`.
pub fn any_from_si(data: f64, unit_name: &str) -> Option<f64> {
    from_si(data, unit_name)
}

/// Converts from LS to any unit.
/// OCCT: `UnitsAPI::AnyFromLS(data, unit_name)`.
pub fn any_from_ls(data: f64, unit_name: &str) -> Option<f64> {
    let ls_unit = default_unit_for_quantity(
        dictionary_of_units().find_unit(unit_name)?.quantity_name.as_str(),
        local_system(),
    );
    let si = to_si(data, &ls_unit)?;
    from_si(si, unit_name)
}

/// Converts between any two units.
/// OCCT: `UnitsAPI::AnyToAny(data, from_unit, to_unit)`.
pub fn any_to_any(data: f64, from_unit: &str, to_unit: &str) -> Option<f64> {
    convert(data, from_unit, to_unit)
}

/// Converts from Current System to any unit.
/// OCCT: `UnitsAPI::CurrentToAny(data, quantity, to_unit)`.
pub fn current_to_any(data: f64, quantity: &str, to_unit: &str) -> Option<f64> {
    let cur = current_unit(quantity);
    let si = to_si(data, &cur)?;
    from_si(si, to_unit)
}

/// Converts from any unit to Current System.
/// OCCT: `UnitsAPI::CurrentFromAny(data, quantity, from_unit)`.
pub fn current_from_any(data: f64, quantity: &str, from_unit: &str) -> Option<f64> {
    let si = to_si(data, from_unit)?;
    let cur = current_unit(quantity);
    from_si(si, &cur)
}

/// OCCT: `UnitsAPI::Check(quantity, unit)` — verify a unit belongs to a quantity.
pub fn check(quantity: &str, unit_name: &str) -> bool {
    let dict = dictionary_of_units();
    if let Some(u) = dict.find_unit(unit_name) {
        u.quantity_name.eq_ignore_ascii_case(quantity)
    } else {
        false
    }
}

// ============================================================================
// Convenience for STEP/tooling
// ============================================================================

/// Get the current length unit as a string (e.g. "mm", "m", "in").
/// Uses Current System if set, otherwise LS default.
pub fn current_length_unit() -> String {
    current_unit("Length")
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dimensions_length() {
        let d = Dimensions::length();
        assert!((d.length - 1.0).abs() < 1e-12);
        assert!((d.mass - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dimensions_mul_div_pow() {
        let l = Dimensions::length();
        let t = Dimensions::time();
        let vel = l.multiply(&t.power(-1.0));
        // Velocity = L * T^(-1) → length=1, time=-1
        assert!((vel.length - 1.0).abs() < 1e-12);
        assert!((vel.time - (-1.0)).abs() < 1e-12);

        let area = l.multiply(&l);
        assert!((area.length - 2.0).abs() < 1e-12);

        let back = vel.divide(&t.power(-1.0));
        assert!((back.length - 1.0).abs() < 1e-12);
        assert!((back.time - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_dimensions_is_equal() {
        let a = Dimensions::length();
        let b = Dimensions::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!(a.is_equal(&b));
        assert!(!a.is_equal(&Dimensions::mass()));
    }

    #[test]
    fn test_dictionary_built() {
        let dict = dictionary_of_units();
        assert!(dict.quantities().len() >= 9);
        assert!(dict.find_quantity("Length").is_some());
        assert!(dict.find_quantity("Mass").is_some());
        assert!(dict.find_quantity("PlaneAngle").is_some());
    }

    #[test]
    fn test_find_unit() {
        let dict = dictionary_of_units();
        assert!(dict.find_unit("mm").is_some());
        assert!(dict.find_unit("MM").is_some());
        assert!(dict.find_unit("in").is_some());
        assert!(dict.find_unit("deg").is_some());
        assert!(dict.find_unit("kg").is_some());
        assert!(dict.find_unit("unknown").is_none());
    }

    #[test]
    fn test_to_si() {
        assert!((to_si(1000.0, "mm").unwrap() - 1.0).abs() < 1e-12);
        assert!((to_si(1.0, "in").unwrap() - 0.0254).abs() < 1e-12);
        assert!((to_si(180.0, "deg").unwrap() - std::f64::consts::PI).abs() < 1e-12);
        assert!(to_si(1.0, "unknown").is_none());
    }

    #[test]
    fn test_from_si() {
        assert!((from_si(1.0, "mm").unwrap() - 1000.0).abs() < 1e-12);
        assert!((from_si(0.0254, "in").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_convert() {
        let mm = convert(1.0, "in", "mm").unwrap();
        assert!((mm - 25.4).abs() < 1e-12);
        let inch = convert(25.4, "mm", "in").unwrap();
        assert!((inch - 1.0).abs() < 1e-12);
        // Dimension mismatch
        assert!(convert(1.0, "mm", "deg").is_none());
    }

    #[test]
    fn test_local_system_default() {
        assert_eq!(local_system(), SystemUnits::MDTV);
        set_local_system(SystemUnits::SI);
        assert_eq!(local_system(), SystemUnits::SI);
        set_local_system(SystemUnits::MDTV);
    }

    #[test]
    fn test_ls_to_si() {
        // MDTV default: length unit = mm
        assert!((ls_to_si(1000.0, "Length").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_si_to_ls() {
        // MDTV default: length unit = mm
        assert!((si_to_ls(1.0, "Length").unwrap() - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn test_current_unit() {
        // Default (MDTV): length → mm
        assert_eq!(current_unit("Length"), "mm");
        assert_eq!(current_length_unit(), "mm");

        set_current_unit("Length", "in");
        assert_eq!(current_unit("Length"), "in");

        // Test conversion
        let si = current_to_si(1.0, "Length").unwrap();
        assert!((si - 0.0254).abs() < 1e-12);
        let back = current_from_si(si, "Length").unwrap();
        assert!((back - 1.0).abs() < 1e-12);

        // Reset
        set_current_unit("Length", "mm");
    }

    #[test]
    fn test_any_to_any() {
        let r = any_to_any(1.0, "in", "mm").unwrap();
        assert!((r - 25.4).abs() < 1e-12);
    }

    #[test]
    fn test_check() {
        assert!(check("Length", "mm"));
        assert!(check("Length", "in"));
        assert!(!check("Length", "deg"));
        assert!(!check("Length", "unknown"));
    }

    #[test]
    fn test_all_units_roundtrip() {
        let dict = dictionary_of_units();
        for q in dict.quantities() {
            for u in &q.units {
                for sym in &u.symbols {
                    let si = to_si(1.0, sym).expect(&format!("to_si failed for {}", sym));
                    let back = from_si(si, sym).expect(&format!("from_si failed for {}", sym));
                    assert!(
                        (back - 1.0).abs() < 1e-12,
                        "roundtrip failed for {}: {} → {} → {}",
                        sym,
                        1.0,
                        si,
                        back
                    );
                }
            }
        }
    }
}
