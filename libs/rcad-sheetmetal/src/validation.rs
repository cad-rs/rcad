//! Design for Manufacturing (DFM) validation for sheet metal parts.
//!
//! This module provides validation checks for sheet metal designs to ensure
//! manufacturability and identify potential issues before production.

use glam::DVec2;
use serde::{Deserialize, Serialize};

use crate::features::{SheetMetalPart, SheetMetalMaterial, FlangeFeature, BendFeature, CounterSinkFeature, EmbossFeature, LouverFeature};
use crate::unfold::FlatPattern;

/// Severity level for validation issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Severity {
    /// Critical issue - part cannot be manufactured.
    Critical,
    /// Warning - part may have quality issues.
    Warning,
    /// Information - best practice recommendation.
    Info,
}

impl std::fmt::Display for Severity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Critical => write!(f, "CRITICAL"),
            Self::Warning => write!(f, "WARNING"),
            Self::Info => write!(f, "INFO"),
        }
    }
}

/// Category of validation issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IssueCategory {
    /// Bend-related issues.
    Bend,
    /// Feature spacing issues.
    Spacing,
    /// Dimension issues.
    Dimension,
    /// Material issues.
    Material,
    /// Feature interaction issues.
    FeatureInteraction,
    /// Tooling constraints.
    Tooling,
}

impl std::fmt::Display for IssueCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bend => write!(f, "Bend"),
            Self::Spacing => write!(f, "Spacing"),
            Self::Dimension => write!(f, "Dimension"),
            Self::Material => write!(f, "Material"),
            Self::FeatureInteraction => write!(f, "Feature Interaction"),
            Self::Tooling => write!(f, "Tooling"),
        }
    }
}

/// A single validation issue.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationIssue {
    /// Severity of the issue.
    pub severity: Severity,
    /// Category of the issue.
    pub category: IssueCategory,
    /// Human-readable message.
    pub message: String,
    /// Location of the issue (if applicable).
    pub location: Option<DVec2>,
    /// Feature index involved (if applicable).
    pub feature_index: Option<usize>,
    /// Suggested fix (if applicable).
    pub suggestion: Option<String>,
}

impl ValidationIssue {
    /// Create a new validation issue.
    pub fn new(severity: Severity, category: IssueCategory, message: &str) -> Self {
        Self {
            severity,
            category,
            message: message.to_string(),
            location: None,
            feature_index: None,
            suggestion: None,
        }
    }

    /// Add location information.
    pub fn with_location(mut self, location: DVec2) -> Self {
        self.location = Some(location);
        self
    }

    /// Add feature index.
    pub fn with_feature(mut self, index: usize) -> Self {
        self.feature_index = Some(index);
        self
    }

    /// Add a suggestion for fixing the issue.
    pub fn with_suggestion(mut self, suggestion: &str) -> Self {
        self.suggestion = Some(suggestion.to_string());
        self
    }

    /// Create a critical issue.
    pub fn critical(category: IssueCategory, message: &str) -> Self {
        Self::new(Severity::Critical, category, message)
    }

    /// Create a warning.
    pub fn warning(category: IssueCategory, message: &str) -> Self {
        Self::new(Severity::Warning, category, message)
    }

    /// Create an info message.
    pub fn info(category: IssueCategory, message: &str) -> Self {
        Self::new(Severity::Info, category, message)
    }
}

/// Result of DFM validation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ValidationResult {
    /// All validation issues found.
    pub issues: Vec<ValidationIssue>,
    /// Whether the part passes validation (no critical issues).
    pub passes: bool,
    /// Number of critical issues.
    pub critical_count: usize,
    /// Number of warnings.
    pub warning_count: usize,
    /// Number of info messages.
    pub info_count: usize,
}

impl ValidationResult {
    /// Create a new empty validation result.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an issue to the result.
    pub fn add_issue(&mut self, issue: ValidationIssue) {
        match issue.severity {
            Severity::Critical => self.critical_count += 1,
            Severity::Warning => self.warning_count += 1,
            Severity::Info => self.info_count += 1,
        }
        self.issues.push(issue);
        self.update_passes();
    }

    /// Add multiple issues to the result.
    pub fn add_issues(&mut self, issues: Vec<ValidationIssue>) {
        for issue in issues {
            self.add_issue(issue);
        }
    }

    /// Update the passes flag based on issues.
    fn update_passes(&mut self) {
        self.passes = self.critical_count == 0;
    }

    /// Get all critical issues.
    pub fn critical_issues(&self) -> Vec<&ValidationIssue> {
        self.issues
            .iter()
            .filter(|i| i.severity == Severity::Critical)
            .collect()
    }

    /// Get all warnings.
    pub fn warnings(&self) -> Vec<&ValidationIssue> {
        self.issues
            .iter()
            .filter(|i| i.severity == Severity::Warning)
            .collect()
    }

    /// Merge another validation result into this one.
    pub fn merge(&mut self, other: ValidationResult) {
        self.issues.extend(other.issues);
        self.critical_count += other.critical_count;
        self.warning_count += other.warning_count;
        self.info_count += other.info_count;
        self.update_passes();
    }
}

/// DFM validator for sheet metal parts.
#[derive(Debug, Clone)]
pub struct DfmValidator {
    /// Minimum bend radius as factor of thickness.
    pub min_bend_radius_factor: f64,
    /// Minimum flange length as factor of thickness.
    pub min_flange_length_factor: f64,
    /// Minimum distance between features.
    pub min_feature_distance: f64,
    /// Minimum distance from feature to bend.
    pub min_bend_distance: f64,
    /// Minimum hole diameter as factor of thickness.
    pub min_hole_diameter_factor: f64,
    /// Minimum distance from hole to edge as factor of thickness.
    pub min_hole_edge_distance_factor: f64,
    /// Maximum bend angle in radians.
    pub max_bend_angle: f64,
    /// Minimum emboss height as factor of thickness.
    pub min_emboss_height_factor: f64,
}

impl Default for DfmValidator {
    fn default() -> Self {
        Self {
            min_bend_radius_factor: 1.0,
            min_flange_length_factor: 4.0,
            min_feature_distance: 2.0,
            min_bend_distance: 3.0,
            min_hole_diameter_factor: 1.0,
            min_hole_edge_distance_factor: 1.5,
            max_bend_angle: std::f64::consts::PI * 0.95, // ~171 degrees
            min_emboss_height_factor: 0.1,
        }
    }
}

impl DfmValidator {
    /// Create a new DFM validator with default settings.
    pub fn new() -> Self {
        Self::default()
    }

    /// Validate a sheet metal part.
    pub fn validate(&self, part: &SheetMetalPart) -> ValidationResult {
        let mut result = ValidationResult::new();

        // Check base feature
        if let Some(base) = &part.base {
            self.validate_base(base, &mut result);
        } else {
            result.add_issue(ValidationIssue::critical(
                IssueCategory::Dimension,
                "No base feature defined",
            ));
        }

        // Validate each feature type
        self.validate_flanges(&part.flanges, part.thickness(), &mut result);
        self.validate_bends(&part.bends, part.thickness(), part.material(), &mut result);
        self.validate_countersinks(&part.countersinks, part.thickness(), &mut result);
        self.validate_embosses(&part.embosses, part.thickness(), &mut result);
        self.validate_louvers(&part.louvers, part.thickness(), &mut result);

        // Check feature spacing
        self.validate_feature_spacing(part, &mut result);

        result
    }

    /// Validate the base feature.
    fn validate_base(&self, base: &crate::features::BaseFeature, result: &mut ValidationResult) {
        // Check dimensions
        if base.width <= 0.0 || base.height <= 0.0 {
            result.add_issue(
                ValidationIssue::critical(IssueCategory::Dimension, "Base dimensions must be positive")
                    .with_suggestion("Set width and height to positive values"),
            );
        }

        if base.thickness <= 0.0 {
            result.add_issue(
                ValidationIssue::critical(IssueCategory::Dimension, "Sheet thickness must be positive")
                    .with_suggestion("Set thickness to a positive value"),
            );
        }

        // Check aspect ratio
        let aspect_ratio = base.width / base.height;
        if aspect_ratio > 10.0 || aspect_ratio < 0.1 {
            result.add_issue(
                ValidationIssue::warning(IssueCategory::Dimension, "Extreme aspect ratio may cause handling issues")
                    .with_suggestion("Consider splitting into multiple parts"),
            );
        }
    }

    /// Validate flange features.
    fn validate_flanges(
        &self,
        flanges: &[FlangeFeature],
        thickness: f64,
        result: &mut ValidationResult,
    ) {
        let min_flange_length = thickness * self.min_flange_length_factor;

        for (idx, flange) in flanges.iter().enumerate() {
            // Check flange length
            if flange.length < min_flange_length {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Dimension, "Flange length below minimum")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Increase flange length to at least {:.2} mm",
                            min_flange_length
                        )),
                );
            }

            // Check bend radius
            let min_radius = thickness * self.min_bend_radius_factor;
            if flange.bend_radius < min_radius {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Bend, "Bend radius below minimum")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Increase bend radius to at least {:.2} mm",
                            min_radius
                        )),
                );
            }

            // Check bend angle
            if flange.angle > self.max_bend_angle {
                result.add_issue(
                    ValidationIssue::critical(IssueCategory::Bend, "Bend angle exceeds maximum")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Reduce bend angle to at most {:.1} degrees",
                            self.max_bend_angle.to_degrees()
                        )),
                );
            }
        }
    }

    /// Validate bend features.
    fn validate_bends(
        &self,
        bends: &[BendFeature],
        thickness: f64,
        material: Option<&SheetMetalMaterial>,
        result: &mut ValidationResult,
    ) {
        let min_radius = thickness * self.min_bend_radius_factor;

        for (idx, bend) in bends.iter().enumerate() {
            // Check bend radius
            if bend.radius < min_radius {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Bend, "Bend radius below minimum")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Increase bend radius to at least {:.2} mm",
                            min_radius
                        )),
                );
            }

            // Check for potential material cracking (tight radius with certain materials)
            if let Some(mat) = material {
                if bend.radius < thickness && mat.yield_strength > 300.0e6 {
                    result.add_issue(
                        ValidationIssue::warning(IssueCategory::Material, "Tight bend radius may cause cracking in high-strength material")
                            .with_feature(idx)
                            .with_suggestion(&format!(
                                "Consider using larger bend radius or annealed material"
                            )),
                    );
                }
            }

            // Check bend length
            if bend.length() < thickness {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Dimension, "Very short bend line")
                        .with_feature(idx),
                );
            }
        }
    }

    /// Validate countersink features.
    fn validate_countersinks(
        &self,
        countersinks: &[CounterSinkFeature],
        thickness: f64,
        result: &mut ValidationResult,
    ) {
        let min_hole_diameter = thickness * self.min_hole_diameter_factor;

        for (idx, cs) in countersinks.iter().enumerate() {
            // Check hole diameter
            if cs.hole_diameter < min_hole_diameter {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Tooling, "Hole diameter below minimum")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Increase hole diameter to at least {:.2} mm",
                            min_hole_diameter
                        )),
                );
            }

            // Check countersink depth vs thickness
            let cs_depth = cs.calculated_depth();
            if cs_depth > thickness * 0.8 {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Dimension, "Countersink depth exceeds 80% of thickness")
                        .with_feature(idx)
                        .with_suggestion("Reduce countersink diameter or use thinner material"),
                );
            }

            // Check for countersink on too thin material
            if thickness < 1.0 {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Tooling, "Countersink on thin material may cause distortion")
                        .with_feature(idx),
                );
            }
        }
    }

    /// Validate emboss features.
    fn validate_embosses(
        &self,
        embosses: &[EmbossFeature],
        thickness: f64,
        result: &mut ValidationResult,
    ) {
        let min_emboss_height = thickness * self.min_emboss_height_factor;

        for (idx, emboss) in embosses.iter().enumerate() {
            // Check emboss height
            if emboss.height.abs() < min_emboss_height {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Dimension, "Emboss height very small")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Increase emboss height to at least {:.2} mm",
                            min_emboss_height
                        )),
                );
            }

            // Check emboss aspect ratio
            if emboss.font_size > 0.0 {
                let ratio = emboss.height.abs() / emboss.font_size;
                if ratio > 0.5 {
                    result.add_issue(
                        ValidationIssue::warning(IssueCategory::Dimension, "High emboss aspect ratio may cause tearing")
                            .with_feature(idx)
                            .with_suggestion("Reduce emboss height relative to font size"),
                    );
                }
            }

            // Check for text length
            if emboss.text.len() > 50 {
                result.add_issue(
                    ValidationIssue::info(IssueCategory::Dimension, "Long emboss text may be difficult to read")
                        .with_feature(idx),
                );
            }
        }
    }

    /// Validate louver features.
    fn validate_louvers(
        &self,
        louvers: &[LouverFeature],
        thickness: f64,
        result: &mut ValidationResult,
    ) {
        for (idx, louver) in louvers.iter().enumerate() {
            // Check louver dimensions
            if louver.length < louver.width * 2.0 {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Dimension, "Louver length should be at least 2x width")
                        .with_feature(idx),
                );
            }

            // Check height vs thickness
            if louver.height > thickness * 5.0 {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Dimension, "High louver may cause distortion")
                        .with_feature(idx)
                        .with_suggestion("Reduce louver height or increase material thickness"),
                );
            }

            // Check array spacing
            if louver.count > 1 && louver.spacing < louver.width * 1.5 {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Spacing, "Louver array spacing too small")
                        .with_feature(idx)
                        .with_suggestion(&format!(
                            "Increase spacing to at least {:.2} mm",
                            louver.width * 1.5
                        )),
                );
            }
        }
    }

    /// Validate spacing between features.
    fn validate_feature_spacing(&self, part: &SheetMetalPart, result: &mut ValidationResult) {
        // Check distance between bends and other features
        let thickness = part.thickness();

        // Simple check: warn if too many features in small area
        let total_features = part.feature_count();
        if let Some(base) = &part.base {
            let area = base.area();
            let feature_density = total_features as f64 / area * 10000.0; // per 100mm x 100mm

            if feature_density > 10.0 {
                result.add_issue(
                    ValidationIssue::warning(IssueCategory::Spacing, "High feature density may cause tooling issues")
                        .with_suggestion("Consider spreading features or using multiple operations"),
                );
            }
        }

        // Check countersink proximity to bends
        for (cs_idx, cs) in part.countersinks.iter().enumerate() {
            for (bend_idx, bend) in part.bends.iter().enumerate() {
                let cs_pos = DVec2::new(cs.center.x, cs.center.y);
                let bend_start = DVec2::new(bend.start.x, bend.start.y);

                // Simple distance check
                let dist = (cs_pos - bend_start).length();
                if dist < thickness * self.min_bend_distance {
                    result.add_issue(
                        ValidationIssue::warning(IssueCategory::FeatureInteraction, "Countersink too close to bend")
                            .with_feature(cs_idx)
                            .with_suggestion(&format!(
                                "Move countersink at least {:.2} mm from bend line",
                                thickness * self.min_bend_distance
                            )),
                    );
                }
            }
        }
    }

    /// Validate a flat pattern.
    pub fn validate_flat_pattern(&self, pattern: &FlatPattern) -> ValidationResult {
        let mut result = ValidationResult::new();

        // Check for overlapping bend lines
        for (i, bl1) in pattern.bend_lines.iter().enumerate() {
            for (j, bl2) in pattern.bend_lines.iter().enumerate().skip(i + 1) {
                // Check if bend lines are too close
                let dist = (bl1.start - bl2.start).length().min((bl1.end - bl2.end).length());
                if dist < self.min_bend_distance {
                    result.add_issue(
                        ValidationIssue::warning(IssueCategory::Spacing, "Bend lines too close")
                            .with_suggestion(&format!(
                                "Increase distance between bend lines to at least {:.2} mm",
                                self.min_bend_distance
                            )),
                    );
                }
            }
        }

        // Check overall dimensions
        if pattern.length <= 0.0 || pattern.width <= 0.0 {
            result.add_issue(ValidationIssue::critical(
                IssueCategory::Dimension,
                "Invalid flat pattern dimensions",
            ));
        }

        // Check for very small faces
        for (idx, face) in pattern.faces.iter().enumerate() {
            if face.area() < 1.0 {
                result.add_issue(
                    ValidationIssue::info(IssueCategory::Dimension, "Very small face in flat pattern")
                        .with_feature(idx),
                );
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::features::{BaseFeature, SheetMetalMaterial, FlangeFeature, BendFeature, CounterSinkFeature, EmbossFeature, LouverFeature};

    #[test]
    fn validation_issue_creation() {
        let issue = ValidationIssue::critical(IssueCategory::Bend, "Test issue")
            .with_feature(0)
            .with_suggestion("Fix it");

        assert_eq!(issue.severity, Severity::Critical);
        assert_eq!(issue.category, IssueCategory::Bend);
        assert_eq!(issue.message, "Test issue");
        assert_eq!(issue.feature_index, Some(0));
        assert!(issue.suggestion.is_some());
    }

    #[test]
    fn validation_result_empty() {
        let result = ValidationResult::new();
        assert!(result.passes);
        assert_eq!(result.critical_count, 0);
    }

    #[test]
    fn validation_result_add_issue() {
        let mut result = ValidationResult::new();
        result.add_issue(ValidationIssue::warning(IssueCategory::Dimension, "Test"));

        assert!(result.passes);
        assert_eq!(result.warning_count, 1);
    }

    #[test]
    fn validation_result_add_critical() {
        let mut result = ValidationResult::new();
        result.add_issue(ValidationIssue::critical(IssueCategory::Dimension, "Test"));

        assert!(!result.passes);
        assert_eq!(result.critical_count, 1);
    }

    #[test]
    fn validator_default() {
        let validator = DfmValidator::new();
        assert!((validator.min_bend_radius_factor - 1.0).abs() < 1e-6);
        assert!((validator.min_flange_length_factor - 4.0).abs() < 1e-6);
    }

    #[test]
    fn validate_simple_part() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.passes);
    }

    #[test]
    fn validate_part_without_base() {
        let part = SheetMetalPart::new();
        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(!result.passes);
        assert_eq!(result.critical_count, 1);
    }

    #[test]
    fn validate_small_flange() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // Flange length below minimum (4 * 1.5 = 6.0)
        part.add_flange(FlangeFeature::new(0, 3.0, std::f64::consts::FRAC_PI_2, 2.0));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_tight_bend_radius() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // Bend radius below minimum (1 * 1.5 = 1.5)
        part.add_flange(FlangeFeature::new(0, 20.0, std::f64::consts::FRAC_PI_2, 0.5));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_extreme_bend_angle() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // Bend angle near 180 degrees
        part.add_flange(FlangeFeature::new(0, 20.0, 3.1, 2.0));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.critical_count >= 1);
    }

    #[test]
    fn validate_small_countersink() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        part.add_countersink(CounterSinkFeature::new(
            DVec2::new(10.0, 10.0),
            0.5, // Below minimum (1 * 1.5 = 1.5)
            3.0,
            82.0_f64.to_radians(),
        ));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_deep_countersink() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.0);
        // Countersink depth will exceed 80% of thickness
        part.add_countersink(CounterSinkFeature::new(
            DVec2::new(10.0, 10.0),
            2.0,
            8.0, // Large countersink diameter relative to thickness
            std::f64::consts::FRAC_PI_2,
        ));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        // Should have a warning about depth
        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_emboss_height() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // Very small emboss height
        part.add_emboss(EmbossFeature::new(DVec2::new(20.0, 10.0), "TEST", 0.01, 5.0));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_high_aspect_emboss() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // High aspect ratio emboss
        part.add_emboss(EmbossFeature::new(DVec2::new(20.0, 10.0), "TEST", 10.0, 5.0));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_louver_dimensions() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // Louver with length < 2 * width
        part.add_louver(LouverFeature::new(DVec2::new(20.0, 15.0), 5.0, 10.0, 3.0));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_louver_array_spacing() {
        let mut part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        // Louver array with tight spacing
        part.add_louver(LouverFeature::array(DVec2::new(20.0, 15.0), 30.0, 5.0, 3.0, 5, 5.0));

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        assert!(result.warning_count >= 1);
    }

    #[test]
    fn validate_high_feature_density() {
        let mut part = SheetMetalPart::with_base(20.0, 20.0, 1.0);

        // Add many features in small area
        for i in 0..5 {
            for j in 0..5 {
                part.add_countersink(CounterSinkFeature::new(
                    DVec2::new(2.0 + i as f64 * 3.0, 2.0 + j as f64 * 3.0),
                    1.0,
                    2.0,
                    82.0_f64.to_radians(),
                ));
            }
        }

        let validator = DfmValidator::new();
        let result = validator.validate(&part);

        // Should warn about high feature density
        assert!(result.warning_count >= 1 || result.info_count >= 1);
    }

    #[test]
    fn validate_flat_pattern() {
        let part = SheetMetalPart::with_base(100.0, 50.0, 1.5);
        let pattern = FlatPattern::from_part(&part);

        let validator = DfmValidator::new();
        let result = validator.validate_flat_pattern(&pattern);

        assert!(result.passes);
    }

    #[test]
    fn severity_display() {
        assert_eq!(format!("{}", Severity::Critical), "CRITICAL");
        assert_eq!(format!("{}", Severity::Warning), "WARNING");
        assert_eq!(format!("{}", Severity::Info), "INFO");
    }

    #[test]
    fn category_display() {
        assert_eq!(format!("{}", IssueCategory::Bend), "Bend");
        assert_eq!(format!("{}", IssueCategory::Spacing), "Spacing");
    }

    #[test]
    fn validation_result_merge() {
        let mut result1 = ValidationResult::new();
        result1.add_issue(ValidationIssue::critical(IssueCategory::Bend, "Issue 1"));

        let mut result2 = ValidationResult::new();
        result2.add_issue(ValidationIssue::warning(IssueCategory::Dimension, "Issue 2"));

        result1.merge(result2);

        assert_eq!(result1.critical_count, 1);
        assert_eq!(result1.warning_count, 1);
        assert!(!result1.passes);
    }

    #[test]
    fn validation_result_critical_issues() {
        let mut result = ValidationResult::new();
        result.add_issue(ValidationIssue::critical(IssueCategory::Bend, "Critical"));
        result.add_issue(ValidationIssue::warning(IssueCategory::Dimension, "Warning"));

        let critical = result.critical_issues();
        assert_eq!(critical.len(), 1);
    }

    #[test]
    fn validation_result_warnings() {
        let mut result = ValidationResult::new();
        result.add_issue(ValidationIssue::critical(IssueCategory::Bend, "Critical"));
        result.add_issue(ValidationIssue::warning(IssueCategory::Dimension, "Warning 1"));
        result.add_issue(ValidationIssue::warning(IssueCategory::Dimension, "Warning 2"));

        let warnings = result.warnings();
        assert_eq!(warnings.len(), 2);
    }
}
