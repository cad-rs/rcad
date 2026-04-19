//! Sheet metal unfolding (flat pattern generation).
//!
//! This module provides functionality to unfold a bent sheet metal part
//! into a flat pattern, including bend allowance calculations.

use glam::{DAffine2, DMat2, DVec2, DVec3};
use serde::{Deserialize, Serialize};

use crate::features::{BendFeature, FlangeFeature, SheetMetalPart, SheetMetalMaterial};

/// A bend line in the flat pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BendLine {
    /// Start point of the bend line in flat pattern coordinates.
    pub start: DVec2,
    /// End point of the bend line in flat pattern coordinates.
    pub end: DVec2,
    /// Bend angle in radians.
    pub angle: f64,
    /// Inner bend radius.
    pub radius: f64,
    /// Bend direction (up = true, down = false).
    pub bend_up: bool,
    /// Bend allowance for this bend.
    pub allowance: f64,
    /// Index of the corresponding bend feature.
    pub bend_index: Option<usize>,
}

impl BendLine {
    /// Create a new bend line.
    pub fn new(start: DVec2, end: DVec2, angle: f64, radius: f64, bend_up: bool) -> Self {
        Self {
            start,
            end,
            angle,
            radius,
            bend_up,
            allowance: 0.0,
            bend_index: None,
        }
    }

    /// Get the bend line direction vector.
    pub fn direction(&self) -> DVec2 {
        (self.end - self.start).normalize()
    }

    /// Get the bend line length.
    pub fn length(&self) -> f64 {
        (self.end - self.start).length()
    }

    /// Get the perpendicular direction to the bend line.
    pub fn perpendicular(&self) -> DVec2 {
        let dir = self.direction();
        DVec2::new(-dir.y, dir.x)
    }
}

/// A face in the flat pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlatFace {
    /// Face vertices in flat pattern coordinates.
    pub vertices: Vec<DVec2>,
    /// Original face index in the 3D model.
    pub original_face_index: Option<usize>,
    /// Face type.
    pub face_type: FlatFaceType,
}

/// Types of flat pattern faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FlatFaceType {
    /// Base sheet face.
    Base,
    /// Flange face.
    Flange,
    /// Bend face (developable surface).
    Bend,
}

impl FlatFace {
    /// Create a new flat face.
    pub fn new(vertices: Vec<DVec2>, face_type: FlatFaceType) -> Self {
        Self {
            vertices,
            original_face_index: None,
            face_type,
        }
    }

    /// Calculate the face area.
    pub fn area(&self) -> f64 {
        if self.vertices.len() < 3 {
            return 0.0;
        }

        // Shoelace formula for polygon area
        let n = self.vertices.len();
        let mut sum = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            sum += self.vertices[i].x * self.vertices[j].y;
            sum -= self.vertices[j].x * self.vertices[i].y;
        }
        sum.abs() / 2.0
    }

    /// Calculate the centroid of the face.
    pub fn centroid(&self) -> DVec2 {
        if self.vertices.is_empty() {
            return DVec2::ZERO;
        }

        let sum: DVec2 = self.vertices.iter().sum();
        sum / self.vertices.len() as f64
    }
}

/// A complete flat pattern for a sheet metal part.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct FlatPattern {
    /// All faces in the flat pattern.
    pub faces: Vec<FlatFace>,
    /// All bend lines in the flat pattern.
    pub bend_lines: Vec<BendLine>,
    /// Overall bounding box.
    pub bounding_box: (DVec2, DVec2),
    /// Total flat pattern length (X extent).
    pub length: f64,
    /// Total flat pattern width (Y extent).
    pub width: f64,
    /// Material thickness.
    pub thickness: f64,
    /// Material used for the part.
    pub material: Option<SheetMetalMaterial>,
}

impl FlatPattern {
    /// Create a new empty flat pattern.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a flat pattern from a sheet metal part.
    pub fn from_part(part: &SheetMetalPart) -> Self {
        let mut pattern = Self::new();
        pattern.thickness = part.thickness();
        pattern.material = part.material().cloned();

        if let Some(base) = &part.base {
            // Create the base face
            let base_face = FlatFace::new(
                vec![
                    DVec2::new(0.0, 0.0),
                    DVec2::new(base.width, 0.0),
                    DVec2::new(base.width, base.height),
                    DVec2::new(0.0, base.height),
                ],
                FlatFaceType::Base,
            );
            pattern.faces.push(base_face);

            // Process flanges
            for (idx, flange) in part.flanges.iter().enumerate() {
                let (face, bend_line) = Self::process_flange(
                    flange,
                    idx,
                    base.width,
                    base.height,
                    part.thickness(),
                    part.material(),
                );
                pattern.faces.push(face);
                if let Some(bl) = bend_line {
                    pattern.bend_lines.push(bl);
                }
            }

            // Process bends
            for (idx, bend) in part.bends.iter().enumerate() {
                if let Some(bend_line) = Self::process_bend(bend, idx, part.material()) {
                    pattern.bend_lines.push(bend_line);
                }
            }

            // Calculate bounding box
            pattern.calculate_bounding_box();
        }

        pattern
    }

    /// Process a flange feature into a flat face and bend line.
    fn process_flange(
        flange: &FlangeFeature,
        index: usize,
        base_width: f64,
        base_height: f64,
        thickness: f64,
        material: Option<&SheetMetalMaterial>,
    ) -> (FlatFace, Option<BendLine>) {
        // Simplified: assume flange is on the edge specified
        let k_factor = material.map(|m| m.k_factor).unwrap_or(0.44);
        let bend_allowance = flange.angle * (flange.bend_radius + k_factor * thickness);

        // Calculate flange position based on edge index
        let (start, end, flange_start) = match flange.edge_index % 4 {
            0 => {
                // Bottom edge
                let flange_end = DVec2::new(base_width, -flange.length - bend_allowance);
                (
                    DVec2::new(0.0, 0.0),
                    DVec2::new(base_width, 0.0),
                    DVec2::new(0.0, -flange.length - bend_allowance),
                )
            }
            1 => {
                // Right edge
                let flange_end = DVec2::new(base_width + flange.length + bend_allowance, base_height);
                (
                    DVec2::new(base_width, 0.0),
                    DVec2::new(base_width, base_height),
                    DVec2::new(base_width + flange.length + bend_allowance, 0.0),
                )
            }
            2 => {
                // Top edge
                let flange_end = DVec2::new(base_width, base_height + flange.length + bend_allowance);
                (
                    DVec2::new(0.0, base_height),
                    DVec2::new(base_width, base_height),
                    DVec2::new(0.0, base_height + flange.length + bend_allowance),
                )
            }
            _ => {
                // Left edge
                let flange_end = DVec2::new(-flange.length - bend_allowance, base_height);
                (
                    DVec2::new(0.0, 0.0),
                    DVec2::new(0.0, base_height),
                    DVec2::new(-flange.length - bend_allowance, 0.0),
                )
            }
        };

        let face = FlatFace::new(
            vec![start, end, flange_start + (end - start), flange_start],
            FlatFaceType::Flange,
        );

        let bend_line = BendLine::new(start, end, flange.angle, flange.bend_radius, true);
        let mut bend_line = bend_line;
        bend_line.allowance = bend_allowance;
        bend_line.bend_index = Some(index);
        (face, Some(bend_line))
    }

    /// Process a bend feature into a bend line.
    fn process_bend(bend: &BendFeature, index: usize, material: Option<&SheetMetalMaterial>) -> Option<BendLine> {
        let k_factor = material.map(|m| m.k_factor).unwrap_or(0.44);
        let thickness = material.map(|_| 1.0).unwrap_or(1.0); // Default thickness if not available

        let start = DVec2::new(bend.start.x, bend.start.y);
        let end = DVec2::new(bend.end.x, bend.end.y);

        let mut bend_line = BendLine::new(start, end, bend.angle, bend.radius, bend.bend_up);
        bend_line.allowance = bend.angle * (bend.radius + k_factor * thickness);
        bend_line.bend_index = Some(index);

        Some(bend_line)
    }

    /// Calculate the bounding box of the flat pattern.
    fn calculate_bounding_box(&mut self) {
        let mut min = DVec2::new(f64::MAX, f64::MAX);
        let mut max = DVec2::new(f64::MIN, f64::MIN);

        for face in &self.faces {
            for v in &face.vertices {
                min = min.min(*v);
                max = max.max(*v);
            }
        }

        self.bounding_box = (min, max);
        self.length = max.x - min.x;
        self.width = max.y - min.y;
    }

    /// Calculate total material area.
    pub fn total_area(&self) -> f64 {
        self.faces.iter().map(|f| f.area()).sum()
    }

    /// Calculate total material volume.
    pub fn total_volume(&self) -> f64 {
        self.total_area() * self.thickness
    }

    /// Calculate total material mass.
    pub fn total_mass(&self) -> f64 {
        let density = self.material.as_ref().map(|m| m.density).unwrap_or(7850.0);
        self.total_volume() * density
    }

    /// Get bend line at a specific position.
    pub fn get_bend_line_at(&self, position: DVec2, tolerance: f64) -> Option<&BendLine> {
        for bend_line in &self.bend_lines {
            // Check if position is near the bend line
            let dist = point_to_line_distance(position, bend_line.start, bend_line.end);
            if dist <= tolerance {
                return Some(bend_line);
            }
        }
        None
    }

    /// Get all bend lines in a region.
    pub fn get_bend_lines_in_region(&self, min: DVec2, max: DVec2) -> Vec<&BendLine> {
        self.bend_lines
            .iter()
            .filter(|bl| {
                let start_in = bl.start.x >= min.x && bl.start.x <= max.x
                    && bl.start.y >= min.y && bl.start.y <= max.y;
                let end_in = bl.end.x >= min.x && bl.end.x <= max.x
                    && bl.end.y >= min.y && bl.end.y <= max.y;
                start_in || end_in
            })
            .collect()
    }
}

/// Calculate the distance from a point to a line segment.
fn point_to_line_distance(point: DVec2, line_start: DVec2, line_end: DVec2) -> f64 {
    let line_vec = line_end - line_start;
    let point_vec = point - line_start;
    let line_len = line_vec.length();

    if line_len < 1e-10 {
        return point_vec.length();
    }

    let t = (point_vec.dot(line_vec) / line_len).clamp(0.0, line_len) / line_len;
    let projection = line_start + line_vec * t;
    (point - projection).length()
}

/// L-bend (single bend) calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LBendResult {
    /// Flat length of the first leg.
    pub leg1_length: f64,
    /// Flat length of the second leg.
    pub leg2_length: f64,
    /// Bend allowance.
    pub bend_allowance: f64,
    /// Total flat length.
    pub total_flat_length: f64,
}

/// Calculate L-bend flat pattern parameters.
pub fn calculate_l_bend(
    leg1_length: f64,
    leg2_length: f64,
    bend_angle: f64,
    bend_radius: f64,
    thickness: f64,
    k_factor: f64,
) -> LBendResult {
    let bend_allowance = bend_angle * (bend_radius + k_factor * thickness);
    let total_flat_length = leg1_length + leg2_length + bend_allowance;

    LBendResult {
        leg1_length,
        leg2_length,
        bend_allowance,
        total_flat_length,
    }
}

/// U-bend (double bend) calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UBendResult {
    /// Flat length of the first leg.
    pub leg1_length: f64,
    /// Flat length of the middle section.
    pub middle_length: f64,
    /// Flat length of the second leg.
    pub leg2_length: f64,
    /// Bend allowance for each bend.
    pub bend_allowance: f64,
    /// Total flat length.
    pub total_flat_length: f64,
}

/// Calculate U-bend flat pattern parameters.
pub fn calculate_u_bend(
    leg1_length: f64,
    middle_length: f64,
    leg2_length: f64,
    bend_angle: f64,
    bend_radius: f64,
    thickness: f64,
    k_factor: f64,
) -> UBendResult {
    let bend_allowance = bend_angle * (bend_radius + k_factor * thickness);
    let total_flat_length = leg1_length + middle_length + leg2_length + 2.0 * bend_allowance;

    UBendResult {
        leg1_length,
        middle_length,
        leg2_length,
        bend_allowance,
        total_flat_length,
    }
}

/// Z-bend (offset bend) calculation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZBendResult {
    /// Flat length of the first leg.
    pub leg1_length: f64,
    /// Flat length of the offset section.
    pub offset_length: f64,
    /// Flat length of the second leg.
    pub leg2_length: f64,
    /// Bend allowance for each bend.
    pub bend_allowance: f64,
    /// Total flat length.
    pub total_flat_length: f64,
    /// Minimum offset distance.
    pub min_offset: f64,
}

/// Calculate Z-bend flat pattern parameters.
pub fn calculate_z_bend(
    leg1_length: f64,
    offset_length: f64,
    leg2_length: f64,
    bend_angle: f64,
    bend_radius: f64,
    thickness: f64,
    k_factor: f64,
) -> ZBendResult {
    let bend_allowance = bend_angle * (bend_radius + k_factor * thickness);
    let total_flat_length = leg1_length + offset_length + leg2_length + 2.0 * bend_allowance;

    // Minimum offset to avoid tool interference
    let min_offset = (bend_radius + thickness) * 2.0;

    ZBendResult {
        leg1_length,
        offset_length,
        leg2_length,
        bend_allowance,
        total_flat_length,
        min_offset,
    }
}

/// Unfold error type.
#[derive(Debug)]
pub enum UnfoldError {
    /// Invalid bend configuration.
    InvalidBend(String),
    /// Geometry cannot be unfolded.
    NonDevelopable(String),
    /// Feature not found.
    FeatureNotFound(usize),
}

impl std::fmt::Display for UnfoldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidBend(msg) => write!(f, "Invalid bend: {msg}"),
            Self::NonDevelopable(msg) => write!(f, "Non-developable geometry: {msg}"),
            Self::FeatureNotFound(id) => write!(f, "Feature not found: {id}"),
        }
    }
}

impl std::error::Error for UnfoldError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::features::{BaseFeature, FlangeFeature};

    #[test]
    fn bend_line_creation() {
        let bl = BendLine::new(
            DVec2::new(0.0, 0.0),
            DVec2::new(50.0, 0.0),
            std::f64::consts::FRAC_PI_2,
            2.0,
            true,
        );
        assert!((bl.length() - 50.0).abs() < 1e-6);
        assert!((bl.direction().x - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bend_line_perpendicular() {
        let bl = BendLine::new(
            DVec2::new(0.0, 0.0),
            DVec2::new(50.0, 0.0),
            std::f64::consts::FRAC_PI_2,
            2.0,
            true,
        );
        let perp = bl.perpendicular();
        assert!((perp.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn flat_face_area() {
        let face = FlatFace::new(
            vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(10.0, 0.0),
                DVec2::new(10.0, 20.0),
                DVec2::new(0.0, 20.0),
            ],
            FlatFaceType::Base,
        );
        assert!((face.area() - 200.0).abs() < 1e-6);
    }

    #[test]
    fn flat_face_centroid() {
        let face = FlatFace::new(
            vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(10.0, 0.0),
                DVec2::new(10.0, 10.0),
                DVec2::new(0.0, 10.0),
            ],
            FlatFaceType::Base,
        );
        let centroid = face.centroid();
        assert!((centroid.x - 5.0).abs() < 1e-6);
        assert!((centroid.y - 5.0).abs() < 1e-6);
    }

    #[test]
    fn flat_pattern_from_part() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        let pattern = FlatPattern::from_part(&part);
        assert_eq!(pattern.faces.len(), 1);
        assert_eq!(pattern.thickness, 1.5);
    }

    #[test]
    fn flat_pattern_with_flange() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        part.add_flange(FlangeFeature::new(0, 20.0, std::f64::consts::FRAC_PI_2, 2.0));

        let pattern = FlatPattern::from_part(&part);
        assert!(pattern.faces.len() >= 1);
        assert!(!pattern.bend_lines.is_empty());
    }

    #[test]
    fn flat_pattern_bounding_box() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        let pattern = FlatPattern::from_part(&part);

        let (min, max) = pattern.bounding_box;
        assert!((min.x - 0.0).abs() < 1e-6);
        assert!((min.y - 0.0).abs() < 1e-6);
        assert!((max.x - 100.0).abs() < 1e-6);
        assert!((max.y - 50.0).abs() < 1e-6);
    }

    #[test]
    fn flat_pattern_total_area() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        let pattern = FlatPattern::from_part(&part);
        assert!((pattern.total_area() - 5000.0).abs() < 1e-6);
    }

    #[test]
    fn flat_pattern_total_volume() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        let pattern = FlatPattern::from_part(&part);
        assert!((pattern.total_volume() - 7500.0).abs() < 1e-6);
    }

    #[test]
    fn l_bend_calculation() {
        let result = calculate_l_bend(
            50.0,
            30.0,
            std::f64::consts::FRAC_PI_2,
            2.0,
            1.0,
            0.44,
        );
        assert!(result.bend_allowance > 0.0);
        assert!((result.total_flat_length - (50.0 + 30.0 + result.bend_allowance)).abs() < 1e-6);
    }

    #[test]
    fn u_bend_calculation() {
        let result = calculate_u_bend(
            20.0,
            40.0,
            20.0,
            std::f64::consts::FRAC_PI_2,
            2.0,
            1.0,
            0.44,
        );
        assert!(result.bend_allowance > 0.0);
        let expected = 20.0 + 40.0 + 20.0 + 2.0 * result.bend_allowance;
        assert!((result.total_flat_length - expected).abs() < 1e-6);
    }

    #[test]
    fn z_bend_calculation() {
        let result = calculate_z_bend(
            30.0,
            10.0,
            30.0,
            std::f64::consts::FRAC_PI_2,
            2.0,
            1.0,
            0.44,
        );
        assert!(result.bend_allowance > 0.0);
        assert!(result.min_offset > 0.0);
    }

    #[test]
    fn point_to_line_distance_on_line() {
        let dist = point_to_line_distance(
            DVec2::new(5.0, 0.0),
            DVec2::new(0.0, 0.0),
            DVec2::new(10.0, 0.0),
        );
        assert!((dist - 0.0).abs() < 1e-6);
    }

    #[test]
    fn point_to_line_distance_off_line() {
        let dist = point_to_line_distance(
            DVec2::new(5.0, 5.0),
            DVec2::new(0.0, 0.0),
            DVec2::new(10.0, 0.0),
        );
        assert!((dist - 5.0).abs() < 1e-6);
    }

    #[test]
    fn point_to_line_distance_beyond_segment() {
        let dist = point_to_line_distance(
            DVec2::new(15.0, 5.0),
            DVec2::new(0.0, 0.0),
            DVec2::new(10.0, 0.0),
        );
        // Should be distance to the endpoint
        let expected = (DVec2::new(15.0, 5.0) - DVec2::new(10.0, 0.0)).length();
        assert!((dist - expected).abs() < 1e-6);
    }

    #[test]
    fn get_bend_line_at_position() {
        let mut pattern = FlatPattern::new();
        pattern.bend_lines.push(BendLine::new(
            DVec2::new(0.0, 0.0),
            DVec2::new(50.0, 0.0),
            std::f64::consts::FRAC_PI_2,
            2.0,
            true,
        ));

        let found = pattern.get_bend_line_at(DVec2::new(25.0, 0.5), 1.0);
        assert!(found.is_some());

        let not_found = pattern.get_bend_line_at(DVec2::new(25.0, 5.0), 1.0);
        assert!(not_found.is_none());
    }

    #[test]
    fn get_bend_lines_in_region() {
        let mut pattern = FlatPattern::new();
        pattern.bend_lines.push(BendLine::new(
            DVec2::new(0.0, 0.0),
            DVec2::new(50.0, 0.0),
            std::f64::consts::FRAC_PI_2,
            2.0,
            true,
        ));
        pattern.bend_lines.push(BendLine::new(
            DVec2::new(100.0, 100.0),
            DVec2::new(150.0, 100.0),
            std::f64::consts::FRAC_PI_2,
            2.0,
            true,
        ));

        let in_region = pattern.get_bend_lines_in_region(
            DVec2::new(-10.0, -10.0),
            DVec2::new(60.0, 10.0),
        );
        assert_eq!(in_region.len(), 1);
    }

    #[test]
    fn flat_pattern_mass_calculation() {
        let material = SheetMetalMaterial::new("Steel", 7850.0, 200.0e9, 0.3, 250.0e6, 0.44);
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        if let Some(ref mut base) = part.base {
            base.material = material;
        }

        let pattern = FlatPattern::from_part(&part);
        // Volume = 100 * 50 * 1.5 = 7500 mm^3 = 7.5e-6 m^3
        // Mass = 7.5e-6 * 7850 = 0.058875 kg
        let expected_mass = 7500.0 * 7850.0;
        assert!((pattern.total_mass() - expected_mass).abs() < 1e-3);
    }
}
