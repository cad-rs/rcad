//! Sheet metal feature definitions and operations.
//!
//! This module provides feature primitives for sheet metal design including:
//! - Base feature (starting flat sheet)
//! - Flange feature (bent edges)
//! - Bend feature (generic bends)
//! - Cut feature (cutouts)
//! - Hem feature (edge hems)
//! - CornerRelief feature (corner stress relief)
//! - Louver feature (ventilation louvers)
//! - Dimple feature (dimple formations)
//! - Jog feature (offset bends)
//! - Tab feature (attachment tabs)
//! - CounterSink feature (countersink holes for fasteners)
//! - Emboss feature (embossed text/logos)

use glam::{DAffine3, DMat3, DVec2, DVec3};
use serde::{Deserialize, Serialize};

/// Sheet metal material properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SheetMetalMaterial {
    /// Material name (e.g., "Steel", "Aluminum", "Copper")
    pub name: String,
    /// Material density in kg/m^3
    pub density: f64,
    /// Young's modulus in Pa
    pub youngs_modulus: f64,
    /// Poisson's ratio
    pub poisson_ratio: f64,
    /// Yield strength in Pa
    pub yield_strength: f64,
    /// K-factor for bend allowance calculation (typically 0.3-0.5)
    pub k_factor: f64,
}

impl Default for SheetMetalMaterial {
    fn default() -> Self {
        Self {
            name: "Steel".to_string(),
            density: 7850.0,
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            yield_strength: 250.0e6,
            k_factor: 0.44,
        }
    }
}

impl SheetMetalMaterial {
    /// Create a new sheet metal material with the given properties.
    pub fn new(name: &str, density: f64, youngs_modulus: f64, poisson_ratio: f64, yield_strength: f64, k_factor: f64) -> Self {
        Self {
            name: name.to_string(),
            density,
            youngs_modulus,
            poisson_ratio,
            yield_strength,
            k_factor,
        }
    }

    /// Calculate bend allowance for a given bend angle, radius, and thickness.
    /// Bend allowance = angle * (radius + k_factor * thickness)
    pub fn bend_allowance(&self, angle_rad: f64, inner_radius: f64, thickness: f64) -> f64 {
        angle_rad * (inner_radius + self.k_factor * thickness)
    }

    /// Calculate bend deduction for a given bend.
    /// Bend deduction = 2 * (radius + thickness) * tan(angle/2) - bend_allowance
    pub fn bend_deduction(&self, angle_rad: f64, inner_radius: f64, thickness: f64) -> f64 {
        let ba = self.bend_allowance(angle_rad, inner_radius, thickness);
        2.0 * (inner_radius + thickness) * (angle_rad / 2.0).tan() - ba
    }
}

/// Feature identifier for sheet metal operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct FeatureId(pub usize);

/// Base feature - the starting flat sheet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BaseFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Sheet width (X dimension).
    pub width: f64,
    /// Sheet height (Y dimension).
    pub height: f64,
    /// Sheet thickness.
    pub thickness: f64,
    /// Material properties.
    pub material: SheetMetalMaterial,
    /// Origin point in 3D space.
    pub origin: DVec3,
    /// Normal direction of the sheet.
    pub normal: DVec3,
}

impl BaseFeature {
    /// Create a new base feature (flat sheet).
    pub fn new(width: f64, height: f64, thickness: f64) -> Self {
        Self {
            id: FeatureId(0),
            width,
            height,
            thickness,
            material: SheetMetalMaterial::default(),
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        }
    }

    /// Create with custom material.
    pub fn with_material(mut self, material: SheetMetalMaterial) -> Self {
        self.material = material;
        self
    }

    /// Calculate the flat pattern area.
    pub fn area(&self) -> f64 {
        self.width * self.height
    }

    /// Calculate the volume.
    pub fn volume(&self) -> f64 {
        self.area() * self.thickness
    }

    /// Calculate the mass.
    pub fn mass(&self) -> f64 {
        self.volume() * self.material.density
    }
}

/// Flange feature - a bent edge extending from the base.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlangeFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Edge index on parent feature where flange is attached.
    pub edge_index: usize,
    /// Flange length (extension from bend line).
    pub length: f64,
    /// Bend angle in radians (pi/2 = 90 degrees).
    pub angle: f64,
    /// Inner bend radius.
    pub bend_radius: f64,
    /// Flange offset from edge (0 = flush).
    pub offset: f64,
}

impl FlangeFeature {
    /// Create a new flange feature.
    pub fn new(edge_index: usize, length: f64, angle: f64, bend_radius: f64) -> Self {
        Self {
            id: FeatureId(0),
            edge_index,
            length,
            angle,
            bend_radius,
            offset: 0.0,
        }
    }

    /// Calculate the flange extension length including bend.
    pub fn total_length(&self, thickness: f64, k_factor: f64) -> f64 {
        let bend_allowance = self.angle * (self.bend_radius + k_factor * thickness);
        self.length + bend_allowance
    }
}

/// Bend feature - a generic bend in the sheet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BendFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Start point of bend line.
    pub start: DVec3,
    /// End point of bend line.
    pub end: DVec3,
    /// Bend angle in radians.
    pub angle: f64,
    /// Inner bend radius.
    pub radius: f64,
    /// Direction of bend (up or down relative to sheet normal).
    pub bend_up: bool,
}

impl BendFeature {
    /// Create a new bend feature.
    pub fn new(start: DVec3, end: DVec3, angle: f64, radius: f64) -> Self {
        Self {
            id: FeatureId(0),
            start,
            end,
            angle,
            radius,
            bend_up: true,
        }
    }

    /// Get the bend line direction vector.
    pub fn direction(&self) -> DVec3 {
        (self.end - self.start).normalize()
    }

    /// Get the bend line length.
    pub fn length(&self) -> f64 {
        (self.end - self.start).length()
    }
}

/// Cut feature - a cutout in the sheet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CutFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Cut profile vertices (2D coordinates on sheet).
    pub profile: Vec<DVec2>,
    /// Cut depth (0 = through cut).
    pub depth: f64,
    /// Is this an internal cut (true) or external profile cut (false).
    pub internal: bool,
}

impl CutFeature {
    /// Create a new through cut.
    pub fn new(profile: Vec<DVec2>) -> Self {
        Self {
            id: FeatureId(0),
            profile,
            depth: 0.0,
            internal: true,
        }
    }

    /// Create a blind cut (partial depth).
    pub fn blind(profile: Vec<DVec2>, depth: f64) -> Self {
        Self {
            id: FeatureId(0),
            profile,
            depth,
            internal: true,
        }
    }
}

/// Hem feature - a folded edge for safety or edge finishing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HemFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Edge index where hem is applied.
    pub edge_index: usize,
    /// Hem type.
    pub hem_type: HemType,
    /// Hem width.
    pub width: f64,
}

/// Types of hem features.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HemType {
    /// Open hem (gap between folded portion and sheet).
    Open,
    /// Closed hem (folded portion touches sheet).
    Closed,
    /// Teardrop hem (curved fold).
    Teardrop,
}

impl HemFeature {
    /// Create a new hem feature.
    pub fn new(edge_index: usize, hem_type: HemType, width: f64) -> Self {
        Self {
            id: FeatureId(0),
            edge_index,
            hem_type,
            width,
        }
    }
}

/// Corner relief feature - stress relief at bent corners.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CornerReliefFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Corner position.
    pub corner: DVec2,
    /// Relief type.
    pub relief_type: CornerReliefType,
    /// Relief size (radius or width).
    pub size: f64,
}

/// Types of corner relief.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CornerReliefType {
    /// Circular relief.
    Circular,
    /// Square relief.
    Square,
    /// V-notch relief.
    VNotch,
}

impl CornerReliefFeature {
    /// Create a new corner relief feature.
    pub fn new(corner: DVec2, relief_type: CornerReliefType, size: f64) -> Self {
        Self {
            id: FeatureId(0),
            corner,
            relief_type,
            size,
        }
    }
}

/// Louver feature - ventilation louvers for airflow.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LouverFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Center position of the louver.
    pub center: DVec2,
    /// Louver length.
    pub length: f64,
    /// Louver width.
    pub width: f64,
    /// Louver height (how much it protrudes).
    pub height: f64,
    /// Opening angle in radians.
    pub angle: f64,
    /// Number of louvers in the array (1 for single).
    pub count: usize,
    /// Spacing between louvers (for arrays).
    pub spacing: f64,
}

impl LouverFeature {
    /// Create a new single louver.
    pub fn new(center: DVec2, length: f64, width: f64, height: f64) -> Self {
        Self {
            id: FeatureId(0),
            center,
            length,
            width,
            height,
            angle: std::f64::consts::FRAC_PI_4, // 45 degrees default
            count: 1,
            spacing: 0.0,
        }
    }

    /// Create a louver array.
    pub fn array(center: DVec2, length: f64, width: f64, height: f64, count: usize, spacing: f64) -> Self {
        Self {
            id: FeatureId(0),
            center,
            length,
            width,
            height,
            angle: std::f64::consts::FRAC_PI_4,
            count,
            spacing,
        }
    }

    /// Calculate the open area of the louver.
    pub fn open_area(&self) -> f64 {
        let single_area = self.length * self.width * self.angle.sin();
        single_area * self.count as f64
    }
}

/// Dimple feature - small protrusions for stiffening or alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DimpleFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Center position of the dimple.
    pub center: DVec2,
    /// Dimple diameter.
    pub diameter: f64,
    /// Dimple depth/height.
    pub depth: f64,
}

impl DimpleFeature {
    /// Create a new dimple feature.
    pub fn new(center: DVec2, diameter: f64, depth: f64) -> Self {
        Self {
            id: FeatureId(0),
            center,
            diameter,
            depth,
        }
    }
}

/// Jog feature - an offset bend creating a step in the sheet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JogFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Jog line start point.
    pub start: DVec3,
    /// Jog line end point.
    pub end: DVec3,
    /// Jog height (step height).
    pub height: f64,
    /// Jog length (flat section).
    pub length: f64,
    /// Bend radius.
    pub radius: f64,
}

impl JogFeature {
    /// Create a new jog feature.
    pub fn new(start: DVec3, end: DVec3, height: f64, length: f64, radius: f64) -> Self {
        Self {
            id: FeatureId(0),
            start,
            end,
            height,
            length,
            radius,
        }
    }

    /// Get the jog direction.
    pub fn direction(&self) -> DVec3 {
        (self.end - self.start).normalize()
    }
}

/// Tab feature - attachment tabs for fastening.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TabFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Tab position (center of tab base).
    pub position: DVec2,
    /// Tab width.
    pub width: f64,
    /// Tab length (extension from edge).
    pub length: f64,
    /// Tab type.
    pub tab_type: TabType,
}

/// Types of tabs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TabType {
    /// Rectangular tab.
    Rectangular,
    /// Trapezoidal tab (tapered).
    Trapezoidal,
    /// Round tab (semi-circular end).
    Round,
}

impl TabFeature {
    /// Create a new tab feature.
    pub fn new(position: DVec2, width: f64, length: f64, tab_type: TabType) -> Self {
        Self {
            id: FeatureId(0),
            position,
            width,
            length,
            tab_type,
        }
    }
}

/// Counter sink feature - countersink holes for fasteners.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CounterSinkFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Center position of the countersink.
    pub center: DVec2,
    /// Main hole diameter.
    pub hole_diameter: f64,
    /// Countersink diameter.
    pub countersink_diameter: f64,
    /// Countersink angle in radians (typically 60, 82, 90, 100, 120, or 150 degrees).
    pub angle: f64,
    /// Countersink depth (calculated if 0).
    pub depth: f64,
}

impl CounterSinkFeature {
    /// Create a new countersink feature.
    pub fn new(center: DVec2, hole_diameter: f64, countersink_diameter: f64, angle: f64) -> Self {
        Self {
            id: FeatureId(0),
            center,
            hole_diameter,
            countersink_diameter,
            angle,
            depth: 0.0,
        }
    }

    /// Create a countersink with a specific depth.
    pub fn with_depth(center: DVec2, hole_diameter: f64, countersink_diameter: f64, angle: f64, depth: f64) -> Self {
        Self {
            id: FeatureId(0),
            center,
            hole_diameter,
            countersink_diameter,
            angle,
            depth,
        }
    }

    /// Calculate the countersink depth if not specified.
    /// Depth = (countersink_diameter - hole_diameter) / (2 * tan(angle/2))
    pub fn calculated_depth(&self) -> f64 {
        if self.depth > 0.0 {
            self.depth
        } else {
            (self.countersink_diameter - self.hole_diameter) / (2.0 * (self.angle / 2.0).tan())
        }
    }

    /// Calculate the area of material removed.
    pub fn removed_area(&self) -> f64 {
        let hole_area = std::f64::consts::PI * (self.hole_diameter / 2.0).powi(2);
        let countersink_area = std::f64::consts::PI * (self.countersink_diameter / 2.0).powi(2);
        countersink_area - hole_area
    }

    /// Create a standard 82-degree countersink (most common).
    pub fn standard_82_degree(center: DVec2, hole_diameter: f64, countersink_diameter: f64) -> Self {
        Self::new(center, hole_diameter, countersink_diameter, 82.0_f64.to_radians())
    }

    /// Create a standard 90-degree countersink.
    pub fn standard_90_degree(center: DVec2, hole_diameter: f64, countersink_diameter: f64) -> Self {
        Self::new(center, hole_diameter, countersink_diameter, std::f64::consts::FRAC_PI_2)
    }
}

/// Emboss feature - embossed text or logos on the sheet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmbossFeature {
    /// Unique identifier for this feature.
    pub id: FeatureId,
    /// Center position of the emboss area.
    pub center: DVec2,
    /// Text or identifier for the emboss pattern.
    pub text: String,
    /// Emboss height (positive = raised, negative = debossed).
    pub height: f64,
    /// Font size (approximate).
    pub font_size: f64,
    /// Width of the emboss area.
    pub width: f64,
    /// Height of the emboss area (bounding box).
    pub bounding_height: f64,
    /// Font style.
    pub font_style: FontStyle,
}

/// Font styles for emboss features.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FontStyle {
    /// Regular/plain style.
    Regular,
    /// Bold style.
    Bold,
    /// Italic style.
    Italic,
}

impl EmbossFeature {
    /// Create a new text emboss feature.
    pub fn new(center: DVec2, text: &str, height: f64, font_size: f64) -> Self {
        let estimated_width = text.len() as f64 * font_size * 0.6;
        Self {
            id: FeatureId(0),
            center,
            text: text.to_string(),
            height,
            font_size,
            width: estimated_width,
            bounding_height: font_size,
            font_style: FontStyle::Regular,
        }
    }

    /// Create a deboss (recessed) feature.
    pub fn deboss(center: DVec2, text: &str, depth: f64, font_size: f64) -> Self {
        let mut feature = Self::new(center, text, -depth.abs(), font_size);
        feature.width = text.len() as f64 * font_size * 0.6;
        feature
    }

    /// Set the font style.
    pub fn with_style(mut self, style: FontStyle) -> Self {
        self.font_style = style;
        self
    }

    /// Calculate the bounding box of the emboss area.
    pub fn bounding_box(&self) -> (DVec2, DVec2) {
        let half_w = self.width / 2.0;
        let half_h = self.bounding_height / 2.0;
        (
            DVec2::new(self.center.x - half_w, self.center.y - half_h),
            DVec2::new(self.center.x + half_w, self.center.y + half_h),
        )
    }

    /// Check if this is a raised emboss.
    pub fn is_raised(&self) -> bool {
        self.height > 0.0
    }

    /// Check if this is a deboss (recessed).
    pub fn is_deboss(&self) -> bool {
        self.height < 0.0
    }
}

/// Error type for sheet metal feature operations.
#[derive(Debug)]
pub enum FeatureError {
    /// Invalid dimensions provided.
    InvalidDimensions(String),
    /// Feature not found.
    FeatureNotFound(FeatureId),
    /// Invalid geometry.
    InvalidGeometry(String),
    /// Feature conflict.
    FeatureConflict(String),
}

impl std::fmt::Display for FeatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDimensions(msg) => write!(f, "Invalid dimensions: {msg}"),
            Self::FeatureNotFound(id) => write!(f, "Feature not found: {:?}", id),
            Self::InvalidGeometry(msg) => write!(f, "Invalid geometry: {msg}"),
            Self::FeatureConflict(msg) => write!(f, "Feature conflict: {msg}"),
        }
    }
}

impl std::error::Error for FeatureError {}

/// Sheet metal part containing all features.
#[derive(Debug, Clone, Default)]
pub struct SheetMetalPart {
    /// Base feature (the sheet).
    pub base: Option<BaseFeature>,
    /// All flange features.
    pub flanges: Vec<FlangeFeature>,
    /// All bend features.
    pub bends: Vec<BendFeature>,
    /// All cut features.
    pub cuts: Vec<CutFeature>,
    /// All hem features.
    pub hems: Vec<HemFeature>,
    /// All corner relief features.
    pub corner_reliefs: Vec<CornerReliefFeature>,
    /// All louver features.
    pub louvers: Vec<LouverFeature>,
    /// All dimple features.
    pub dimples: Vec<DimpleFeature>,
    /// All jog features.
    pub jogs: Vec<JogFeature>,
    /// All tab features.
    pub tabs: Vec<TabFeature>,
    /// All countersink features.
    pub countersinks: Vec<CounterSinkFeature>,
    /// All emboss features.
    pub embosses: Vec<EmbossFeature>,
    /// Feature counter for ID generation.
    next_id: usize,
}

impl SheetMetalPart {
    /// Create a new empty sheet metal part.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a sheet metal part with a base feature.
    pub fn with_base(width: f64, height: f64, thickness: f64) -> Self {
        let mut part = Self::new();
        part.set_base(BaseFeature::new(width, height, thickness));
        part
    }

    /// Set the base feature.
    pub fn set_base(&mut self, base: BaseFeature) {
        self.base = Some(base);
    }

    /// Generate a new unique feature ID.
    fn next_feature_id(&mut self) -> FeatureId {
        let id = self.next_id;
        self.next_id += 1;
        FeatureId(id)
    }

    /// Add a flange feature.
    pub fn add_flange(&mut self, mut flange: FlangeFeature) -> FeatureId {
        let id = self.next_feature_id();
        flange.id = id;
        self.flanges.push(flange);
        id
    }

    /// Add a bend feature.
    pub fn add_bend(&mut self, mut bend: BendFeature) -> FeatureId {
        let id = self.next_feature_id();
        bend.id = id;
        self.bends.push(bend);
        id
    }

    /// Add a cut feature.
    pub fn add_cut(&mut self, mut cut: CutFeature) -> FeatureId {
        let id = self.next_feature_id();
        cut.id = id;
        self.cuts.push(cut);
        id
    }

    /// Add a hem feature.
    pub fn add_hem(&mut self, mut hem: HemFeature) -> FeatureId {
        let id = self.next_feature_id();
        hem.id = id;
        self.hems.push(hem);
        id
    }

    /// Add a corner relief feature.
    pub fn add_corner_relief(&mut self, mut relief: CornerReliefFeature) -> FeatureId {
        let id = self.next_feature_id();
        relief.id = id;
        self.corner_reliefs.push(relief);
        id
    }

    /// Add a louver feature.
    pub fn add_louver(&mut self, mut louver: LouverFeature) -> FeatureId {
        let id = self.next_feature_id();
        louver.id = id;
        self.louvers.push(louver);
        id
    }

    /// Add a dimple feature.
    pub fn add_dimple(&mut self, mut dimple: DimpleFeature) -> FeatureId {
        let id = self.next_feature_id();
        dimple.id = id;
        self.dimples.push(dimple);
        id
    }

    /// Add a jog feature.
    pub fn add_jog(&mut self, mut jog: JogFeature) -> FeatureId {
        let id = self.next_feature_id();
        jog.id = id;
        self.jogs.push(jog);
        id
    }

    /// Add a tab feature.
    pub fn add_tab(&mut self, mut tab: TabFeature) -> FeatureId {
        let id = self.next_feature_id();
        tab.id = id;
        self.tabs.push(tab);
        id
    }

    /// Add a countersink feature.
    pub fn add_countersink(&mut self, mut countersink: CounterSinkFeature) -> FeatureId {
        let id = self.next_feature_id();
        countersink.id = id;
        self.countersinks.push(countersink);
        id
    }

    /// Add an emboss feature.
    pub fn add_emboss(&mut self, mut emboss: EmbossFeature) -> FeatureId {
        let id = self.next_feature_id();
        emboss.id = id;
        self.embosses.push(emboss);
        id
    }

    /// Get the total number of features.
    pub fn feature_count(&self) -> usize {
        self.flanges.len()
            + self.bends.len()
            + self.cuts.len()
            + self.hems.len()
            + self.corner_reliefs.len()
            + self.louvers.len()
            + self.dimples.len()
            + self.jogs.len()
            + self.tabs.len()
            + self.countersinks.len()
            + self.embosses.len()
            + if self.base.is_some() { 1 } else { 0 }
    }

    /// Get the thickness of the base sheet.
    pub fn thickness(&self) -> f64 {
        self.base.as_ref().map(|b| b.thickness).unwrap_or(0.0)
    }

    /// Get the material of the base sheet.
    pub fn material(&self) -> Option<&SheetMetalMaterial> {
        self.base.as_ref().map(|b| &b.material)
    }

    /// Calculate total bend allowance for all bends.
    pub fn total_bend_allowance(&self) -> f64 {
        let k_factor = self.material().map(|m| m.k_factor).unwrap_or(0.44);
        let thickness = self.thickness();

        let from_bends: f64 = self.bends.iter().map(|b| {
            b.angle * (b.radius + k_factor * thickness)
        }).sum();

        let from_flanges: f64 = self.flanges.iter().map(|f| {
            f.angle * (f.bend_radius + k_factor * thickness)
        }).sum();

        from_bends + from_flanges
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_feature_creation() {
        let base = BaseFeature::new(100.0, 50.0, 1.5);
        assert_eq!(base.width, 100.0);
        assert_eq!(base.height, 50.0);
        assert_eq!(base.thickness, 1.5);
        assert!((base.area() - 5000.0).abs() < 1e-6);
        assert!((base.volume() - 7500.0).abs() < 1e-6);
    }

    #[test]
    fn base_feature_with_material() {
        let material = SheetMetalMaterial::new("Aluminum", 2700.0, 70.0e9, 0.33, 280.0e6, 0.42);
        let base = BaseFeature::new(100.0, 100.0, 2.0).with_material(material);
        assert_eq!(base.material.name, "Aluminum");
        assert!((base.mass() - 10000.0 * 2.0 * 2700.0).abs() < 1e-3);
    }

    #[test]
    fn flange_feature_creation() {
        let flange = FlangeFeature::new(0, 25.0, std::f64::consts::FRAC_PI_2, 2.0);
        assert_eq!(flange.edge_index, 0);
        assert_eq!(flange.length, 25.0);
        assert!((flange.angle - std::f64::consts::FRAC_PI_2).abs() < 1e-6);
        assert_eq!(flange.bend_radius, 2.0);
    }

    #[test]
    fn bend_feature_creation() {
        let bend = BendFeature::new(
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(50.0, 0.0, 0.0),
            std::f64::consts::FRAC_PI_4,
            3.0,
        );
        assert!((bend.length() - 50.0).abs() < 1e-6);
        assert!((bend.direction().x - 1.0).abs() < 1e-6);
    }

    #[test]
    fn cut_feature_creation() {
        let cut = CutFeature::new(vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(10.0, 0.0),
            DVec2::new(10.0, 10.0),
            DVec2::new(0.0, 10.0),
        ]);
        assert_eq!(cut.profile.len(), 4);
        assert_eq!(cut.depth, 0.0); // Through cut
        assert!(cut.internal);
    }

    #[test]
    fn blind_cut_feature() {
        let cut = CutFeature::blind(
            vec![DVec2::new(0.0, 0.0), DVec2::new(5.0, 0.0), DVec2::new(2.5, 5.0)],
            1.5,
        );
        assert_eq!(cut.depth, 1.5);
    }

    #[test]
    fn hem_feature_creation() {
        let hem = HemFeature::new(1, HemType::Closed, 5.0);
        assert_eq!(hem.edge_index, 1);
        assert_eq!(hem.hem_type, HemType::Closed);
        assert_eq!(hem.width, 5.0);
    }

    #[test]
    fn corner_relief_creation() {
        let relief = CornerReliefFeature::new(
            DVec2::new(10.0, 10.0),
            CornerReliefType::Circular,
            3.0,
        );
        assert_eq!(relief.size, 3.0);
        assert_eq!(relief.relief_type, CornerReliefType::Circular);
    }

    #[test]
    fn louver_feature_creation() {
        let louver = LouverFeature::new(DVec2::new(20.0, 15.0), 30.0, 5.0, 3.0);
        assert_eq!(louver.length, 30.0);
        assert_eq!(louver.count, 1);
    }

    #[test]
    fn louver_array_creation() {
        let louvers = LouverFeature::array(DVec2::new(0.0, 0.0), 25.0, 4.0, 2.5, 5, 10.0);
        assert_eq!(louvers.count, 5);
        assert_eq!(louvers.spacing, 10.0);
    }

    #[test]
    fn louver_open_area() {
        let louver = LouverFeature::new(DVec2::new(0.0, 0.0), 20.0, 5.0, 2.0);
        let area = louver.open_area();
        assert!(area > 0.0);
        assert!(area < 20.0 * 5.0); // Less than the total louver area
    }

    #[test]
    fn dimple_feature_creation() {
        let dimple = DimpleFeature::new(DVec2::new(15.0, 15.0), 8.0, 2.0);
        assert_eq!(dimple.diameter, 8.0);
        assert_eq!(dimple.depth, 2.0);
    }

    #[test]
    fn jog_feature_creation() {
        let jog = JogFeature::new(
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(40.0, 0.0, 0.0),
            5.0,
            10.0,
            2.0,
        );
        assert_eq!(jog.height, 5.0);
        assert_eq!(jog.length, 10.0);
    }

    #[test]
    fn tab_feature_creation() {
        let tab = TabFeature::new(DVec2::new(10.0, 0.0), 8.0, 5.0, TabType::Rectangular);
        assert_eq!(tab.width, 8.0);
        assert_eq!(tab.length, 5.0);
        assert_eq!(tab.tab_type, TabType::Rectangular);
    }

    #[test]
    fn countersink_feature_creation() {
        let cs = CounterSinkFeature::new(DVec2::new(10.0, 10.0), 5.0, 10.0, 82.0_f64.to_radians());
        assert_eq!(cs.hole_diameter, 5.0);
        assert_eq!(cs.countersink_diameter, 10.0);
    }

    #[test]
    fn countersink_calculated_depth() {
        let cs = CounterSinkFeature::new(DVec2::new(0.0, 0.0), 5.0, 10.0, 90.0_f64.to_radians());
        // depth = (10 - 5) / (2 * tan(45)) = 5 / 2 = 2.5
        let depth = cs.calculated_depth();
        assert!((depth - 2.5).abs() < 1e-6);
    }

    #[test]
    fn countersink_standard_82_degree() {
        let cs = CounterSinkFeature::standard_82_degree(DVec2::new(0.0, 0.0), 4.0, 8.0);
        assert!((cs.angle - 82.0_f64.to_radians()).abs() < 1e-6);
    }

    #[test]
    fn countersink_standard_90_degree() {
        let cs = CounterSinkFeature::standard_90_degree(DVec2::new(0.0, 0.0), 4.0, 8.0);
        assert!((cs.angle - std::f64::consts::FRAC_PI_2).abs() < 1e-6);
    }

    #[test]
    fn countersink_removed_area() {
        let cs = CounterSinkFeature::new(DVec2::new(0.0, 0.0), 4.0, 8.0, std::f64::consts::FRAC_PI_2);
        let area = cs.removed_area();
        // Area = pi * 4^2 - pi * 2^2 = pi * 16 - pi * 4 = pi * 12
        assert!((area - std::f64::consts::PI * 12.0).abs() < 1e-6);
    }

    #[test]
    fn emboss_feature_creation() {
        let emboss = EmbossFeature::new(DVec2::new(20.0, 10.0), "LOGO", 0.5, 5.0);
        assert_eq!(emboss.text, "LOGO");
        assert_eq!(emboss.height, 0.5);
        assert!(emboss.is_raised());
        assert!(!emboss.is_deboss());
    }

    #[test]
    fn deboss_feature_creation() {
        let deboss = EmbossFeature::deboss(DVec2::new(10.0, 5.0), "SERIAL", 0.3, 3.0);
        assert_eq!(deboss.text, "SERIAL");
        assert!(deboss.is_deboss());
        assert!(!deboss.is_raised());
    }

    #[test]
    fn emboss_bounding_box() {
        let emboss = EmbossFeature::new(DVec2::new(0.0, 0.0), "TEST", 0.5, 10.0);
        let (min, max) = emboss.bounding_box();
        assert!(min.x < 0.0);
        assert!(max.x > 0.0);
        assert!(min.y < 0.0);
        assert!(max.y > 0.0);
    }

    #[test]
    fn emboss_with_style() {
        let emboss = EmbossFeature::new(DVec2::new(0.0, 0.0), "TEXT", 0.5, 5.0)
            .with_style(FontStyle::Bold);
        assert_eq!(emboss.font_style, FontStyle::Bold);
    }

    #[test]
    fn sheet_metal_part_creation() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        assert!(part.base.is_some());
        assert_eq!(part.thickness(), 1.5);
    }

    #[test]
    fn sheet_metal_part_add_features() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);

        let fid = part.add_flange(FlangeFeature::new(0, 20.0, std::f64::consts::FRAC_PI_2, 2.0));
        assert_eq!(fid.0, 0);

        let bid = part.add_bend(BendFeature::new(DVec3::ZERO, DVec3::X, std::f64::consts::FRAC_PI_4, 2.0));
        assert_eq!(bid.0, 1);

        let cid = part.add_cut(CutFeature::new(vec![DVec2::ZERO]));
        assert_eq!(cid.0, 2);
    }

    #[test]
    fn sheet_metal_part_feature_count() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);

        part.add_flange(FlangeFeature::new(0, 20.0, std::f64::consts::FRAC_PI_2, 2.0));
        part.add_bend(BendFeature::new(DVec3::ZERO, DVec3::X, std::f64::consts::FRAC_PI_4, 2.0));
        part.add_cut(CutFeature::new(vec![DVec2::ZERO]));
        part.add_hem(HemFeature::new(1, HemType::Open, 5.0));
        part.add_louver(LouverFeature::new(DVec2::ZERO, 20.0, 5.0, 2.0));
        part.add_countersink(CounterSinkFeature::new(DVec2::ZERO, 3.0, 6.0, 82.0_f64.to_radians()));
        part.add_emboss(EmbossFeature::new(DVec2::ZERO, "TEST", 0.5, 3.0));

        // 1 base + 6 added features = 7
        assert_eq!(part.feature_count(), 7);
    }

    #[test]
    fn material_bend_allowance() {
        let material = SheetMetalMaterial::default();
        // Bend allowance = angle * (radius + k_factor * thickness)
        let ba = material.bend_allowance(std::f64::consts::FRAC_PI_2, 2.0, 1.0);
        let expected = std::f64::consts::FRAC_PI_2 * (2.0 + 0.44 * 1.0);
        assert!((ba - expected).abs() < 1e-6);
    }

    #[test]
    fn material_bend_deduction() {
        let material = SheetMetalMaterial::default();
        let bd = material.bend_deduction(std::f64::consts::FRAC_PI_2, 2.0, 1.0);
        // Should be positive for typical bends
        assert!(bd > 0.0);
    }

    #[test]
    fn total_bend_allowance_calculation() {
        let mut part = SheetMetalPart::with_base(50.0, 30.0, 1.0);
        part.add_flange(FlangeFeature::new(0, 20.0, std::f64::consts::FRAC_PI_2, 2.0));
        part.add_bend(BendFeature::new(DVec3::ZERO, DVec3::X, std::f64::consts::FRAC_PI_4, 1.5));

        let tba = part.total_bend_allowance();
        assert!(tba > 0.0);
    }
}
