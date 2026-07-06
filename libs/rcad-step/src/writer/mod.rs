use rcad_kernel::appearance::{Color, StepColor};
use crate::StepGeneralProperty;
use crate::{
    StepDatum, StepDatumSystem, StepDimensionalLocation, StepDimensionalSize,
    StepGeometricTolerance, StepGeometricToleranceWithDatumReference, StepKinematicPair,
    StepPropertyDefinitionRepr,
};
use crate::{
    StepDimensionalTolerance, StepToleranceValue, StepPositionTolerance,
    StepOrientationTolerance, StepFormTolerance, StepRunoutTolerance, StepProfileTolerance,
    StepDatumReferenceFrame, StepDatumTarget, StepToleranceZoneDefinitionEnhanced,
    OrientationToleranceType, FormToleranceType, RunoutToleranceType, ProfileToleranceType,
    DatumTargetType,
};
// View and annotation types
use crate::{
    StepView, StepCameraModelD3, StepViewVolume, ViewVolumeType,
    StepNote, StepAnnotationPlane, StepAnnotationOccurrence,
    StepDimensionCurve, StepTerminatorSymbol, StepDatumFeatureCallout,
    TerminatorType,
};
use rcad_kernel::surface_to_bspline;

/// Selects which STEP application protocol to use when writing a file.
///
/// | Variant | `FILE_SCHEMA` token | Typical use |
/// |---------|---------------------|-------------|
/// | `Ap214` | `AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }` | Legacy automotive/CAD interchange |
/// | `Ap242` | `AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 242 1 1 4 }` | Modern MBD/PMI-aware interchange |
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StepProtocol {
    /// ISO 10303-214 "Automotive Design" (default for backward compatibility).
    #[default]
    Ap214,
    /// ISO 10303-242 "Managed Model Based 3D Engineering".
    Ap242,
}
use self::flat::{BRep, Face};
use rcad_kernel::{BSplineCurve2, Curve2d, Curve3, CurveEval, Surface3, topods};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Write;

pub struct ExportSelection<'a> {
    pub selected_faces: &'a [usize],
    pub selected_edges: &'a [usize],
}

/// STEP header section fields similar to OCCT's OCCSTEP* export controls.
#[derive(Debug, Clone, Default)]
pub struct StepHeader {
    pub description: Option<String>,
    pub implementation_level: Option<String>,
    pub model_name: Option<String>,
    pub time_stamp: Option<String>,
    pub author: Option<String>,
    pub organization: Option<String>,
    pub preprocessor_version: Option<String>,
    pub originating_system: Option<String>,
    pub authorization: Option<String>,
}

/// Unified export options for STEP writing (OCCT-style interchange by default).
#[derive(Debug, Clone)]
pub struct StepWriteOptions {
    pub protocol: StepProtocol,
    pub colors: Option<StepColor>,
    pub properties: Vec<StepGeneralProperty>,
    pub ap242_metadata: Option<StepAp242Metadata>,
    pub header: StepHeader,
    /// Include standalone 1D entities as wireframe overlays in full-model export.
    pub export_standalone_wire_overlay: bool,
}

impl Default for StepWriteOptions {
    fn default() -> Self {
        Self {
            protocol: StepProtocol::default(),
            colors: None,
            properties: Vec::new(),
            ap242_metadata: None,
            header: StepHeader::default(),
            export_standalone_wire_overlay: true,
        }
    }
}

/// Optional AP242 metadata entities to be written alongside geometry.
#[derive(Debug, Clone, Default)]
pub struct StepAp242Metadata {
    pub property_definition_representations: Vec<StepPropertyDefinitionRepr>,
    pub dimensional_locations: Vec<StepDimensionalLocation>,
    pub dimensional_sizes: Vec<StepDimensionalSize>,
    pub geometric_tolerances: Vec<StepGeometricTolerance>,
    pub geometric_tolerances_with_datum_references: Vec<StepGeometricToleranceWithDatumReference>,
    pub datums: Vec<StepDatum>,
    pub datum_systems: Vec<StepDatumSystem>,
    pub kinematic_pairs: Vec<StepKinematicPair>,
    // GDT Extended fields
    pub dimensional_tolerances: Vec<StepDimensionalTolerance>,
    pub tolerance_values: Vec<StepToleranceValue>,
    pub position_tolerances: Vec<StepPositionTolerance>,
    pub orientation_tolerances: Vec<StepOrientationTolerance>,
    pub form_tolerances: Vec<StepFormTolerance>,
    pub runout_tolerances: Vec<StepRunoutTolerance>,
    pub profile_tolerances: Vec<StepProfileTolerance>,
    pub datum_reference_frames: Vec<StepDatumReferenceFrame>,
    pub datum_targets: Vec<StepDatumTarget>,
    pub tolerance_zone_definitions_enhanced: Vec<StepToleranceZoneDefinitionEnhanced>,
    // View and annotation fields
    pub views: Vec<StepView>,
    pub cameras: Vec<StepCameraModelD3>,
    pub view_volumes: Vec<StepViewVolume>,
    pub notes: Vec<StepNote>,
    pub annotation_planes: Vec<StepAnnotationPlane>,
    pub annotation_occurrences: Vec<StepAnnotationOccurrence>,
    pub dimension_curves: Vec<StepDimensionCurve>,
    pub terminator_symbols: Vec<StepTerminatorSymbol>,
    pub datum_feature_callouts: Vec<StepDatumFeatureCallout>,
}

pub struct StepWriter;

/// B-rep STEP export aligns with OCCT-style interchange (`SURFACE_CURVE`/PCURVE, radians, si_metre).
struct FaceExportResult {
    face_ids: Vec<u64>,
    used_triangle_fallback: bool,
}

impl StepWriter {
    /// Export BRep to STEP (OCCT-style interchange: radians, si_metre, surface-scoped curves).
    ///
    /// Uses AP214 protocol by default unless overridden via [`write_string_with_options`].
    pub fn write_string(brep: &topods::BRep, selection: ExportSelection<'_>) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                export_standalone_wire_overlay: true,
                ..Default::default()
            },
        )
    }

    /// Export with a single options object containing protocol, header,
    /// color and AP242 metadata controls.
    pub fn write_string_with_options(
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        options: &StepWriteOptions,
    ) -> String {
        let mut writer = Part21Writer::new_with_protocol_and_header(
            options.protocol,
            options.header.clone(),
            options.export_standalone_wire_overlay,
        );
        writer.write_brep_topods(
            brep,
            selection,
            options.colors.as_ref(),
            &options.properties,
            options.ap242_metadata.as_ref(),
        );
        writer.finish()
    }

    /// Export using the specified STEP application protocol.
    ///
    /// Writes the appropriate `FILE_SCHEMA` and `APPLICATION_PROTOCOL_DEFINITION`
    /// for the chosen protocol.
    pub fn write_string_with_protocol(
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol,
                ..Default::default()
            },
        )
    }

    /// Export with additional generic metadata properties.
    ///
    /// Properties are emitted as `GENERAL_PROPERTY` entities.
    pub fn write_string_with_properties(
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol,
                properties: properties.to_vec(),
                ..Default::default()
            },
        )
    }

    /// Export with generic properties plus AP242 metadata entities.
    pub fn write_string_with_ap242_metadata(
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        metadata: &StepAp242Metadata,
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol,
                properties: properties.to_vec(),
                ap242_metadata: Some(metadata.clone()),
                ..Default::default()
            },
        )
    }

    /// Export with per-face / per-solid color information.
    ///
    /// Colors are written as `STYLED_ITEM` + `PRESENTATION_STYLE_ASSIGNMENT`
    /// + `SURFACE_STYLE_USAGE` + `FILL_AREA_STYLE_COLOUR` + `COLOUR_RGB`.
    pub fn write_string_colored(brep: &topods::BRep, colors: &StepColor) -> String {
        Self::write_string_colored_with_protocol(brep, colors, StepProtocol::Ap214)
    }

    /// Export with per-face / per-solid color information and the specified
    /// STEP application protocol.
    pub fn write_string_colored_with_protocol(
        brep: &topods::BRep,
        colors: &StepColor,
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol,
                colors: Some(colors.clone()),
                ..Default::default()
            },
        )
    }

    /// Stream-based export variant that writes UTF-8 STEP content into any sink.
    pub fn write_to<W: Write>(
        sink: &mut W,
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
    ) -> std::io::Result<()> {
        let step = Self::write_string(brep, selection);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with explicit protocol selection.
    pub fn write_to_with_protocol<W: Write>(
        sink: &mut W,
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        protocol: StepProtocol,
    ) -> std::io::Result<()> {
        let step = Self::write_string_with_protocol(brep, selection, protocol);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with generic metadata properties.
    pub fn write_to_with_properties<W: Write>(
        sink: &mut W,
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        protocol: StepProtocol,
    ) -> std::io::Result<()> {
        let step = Self::write_string_with_properties(brep, selection, properties, protocol);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with AP242 metadata entities.
    pub fn write_to_with_ap242_metadata<W: Write>(
        sink: &mut W,
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        metadata: &StepAp242Metadata,
        protocol: StepProtocol,
    ) -> std::io::Result<()> {
        let step =
            Self::write_string_with_ap242_metadata(brep, selection, properties, metadata, protocol);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with a single options object.
    pub fn write_to_with_options<W: Write>(
        sink: &mut W,
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        options: &StepWriteOptions,
    ) -> std::io::Result<()> {
        let step = Self::write_string_with_options(brep, selection, options);
        sink.write_all(step.as_bytes())
    }
}

struct Part21Writer {
    next_id: u64,
    records: Vec<String>,
    vertex_point_ids: HashMap<usize, u64>,
    surface_ids: HashMap<usize, u64>,
    edge_curve_ids: HashMap<usize, u64>,
    seam_edge_curve_ids: HashMap<usize, u64>,
    edge_geometry_ids: HashMap<usize, u64>,
    /// Set by write_brep_topods; read-only during write.
    topods_brep: *const topods::BRep,
    protocol: StepProtocol,
    header: StepHeader,
    export_standalone_wire_overlay: bool,
    strict_plane_closed_ellipse_done: bool,
}

impl Part21Writer {
    fn new_with_protocol_and_header(
        protocol: StepProtocol,
        header: StepHeader,
        export_standalone_wire_overlay: bool,
    ) -> Self {
        Self {
            next_id: 1,
            records: Vec::new(),
            vertex_point_ids: HashMap::new(),
            surface_ids: HashMap::new(),
            edge_curve_ids: HashMap::new(),
            seam_edge_curve_ids: HashMap::new(),
            edge_geometry_ids: HashMap::new(),
            topods_brep: std::ptr::null(),
            protocol,
            header,
            export_standalone_wire_overlay,
            strict_plane_closed_ellipse_done: false,
        }
    }

    /// Accessor for topods::BRep during writing.
    fn tbrep(&self) -> &topods::BRep {
        unsafe { &*self.topods_brep }
    }

    fn finish(self) -> String {
        let schema_token = match self.protocol {
            StepProtocol::Ap214 => {
                "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }"
            }
            StepProtocol::Ap242 => {
                "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 242 1 1 4 }"
            }
        };
        let mut out = String::new();
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        let file_description = self
            .header
            .description
            .as_deref()
            .unwrap_or("RCAD exported geometry");
        let implementation_level = self
            .header
            .implementation_level
            .as_deref()
            .unwrap_or("2;1");
        out.push_str(&format!(
            "FILE_DESCRIPTION(('{}'),'{}');\n",
            esc_step(file_description),
            esc_step(implementation_level)
        ));

        let model_name = self
            .header
            .model_name
            .as_deref()
            .unwrap_or("rcad_export.step");
        let time_stamp = self
            .header
            .time_stamp
            .as_deref()
            .unwrap_or("2026-04-02T00:00:00");
        let author = self.header.author.as_deref().unwrap_or("");
        let organization = self.header.organization.as_deref().unwrap_or("");
        let preprocessor = self
            .header
            .preprocessor_version
            .as_deref()
            .unwrap_or("RCAD");
        let originating = self
            .header
            .originating_system
            .as_deref()
            .unwrap_or("RCAD");
        let authorization = self.header.authorization.as_deref().unwrap_or("");
        out.push_str(&format!(
            "FILE_NAME('{}','{}',('{}'),('{}'),'{}','{}','{}');\n",
            esc_step(model_name),
            esc_step(time_stamp),
            esc_step(author),
            esc_step(organization),
            esc_step(preprocessor),
            esc_step(originating),
            esc_step(authorization)
        ));
        out.push_str(&format!("FILE_SCHEMA(('{}'));\n", schema_token));
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");
        for record in self.records {
            out.push_str(&record);
            out.push('\n');
        }
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        out
    }

    /// Entry point for topods::BRep — builds FlatBRep internally.
    fn write_brep_topods(
        &mut self,
        brep: &topods::BRep,
        selection: ExportSelection<'_>,
        colors: Option<&StepColor>,
        properties: &[StepGeneralProperty],
        ap242_metadata: Option<&StepAp242Metadata>,
    ) {
        self.topods_brep = brep as *const topods::BRep;
        self.write_brep( selection, colors, properties, ap242_metadata);
    }

    fn write_brep(
        &mut self,
        selection: ExportSelection<'_>,
        colors: Option<&StepColor>,
        properties: &[StepGeneralProperty],
        ap242_metadata: Option<&StepAp242Metadata>,
    ) {
        // Optional general metadata properties.
        for prop in properties {
            let desc = prop.description.as_deref().unwrap_or("");
            let gp = self.general_property(&prop.name, desc);
            self.property_definition(&prop.name, desc, gp);
        }
        if matches!(self.protocol, StepProtocol::Ap242)
            && let Some(meta) = ap242_metadata {
                self.write_ap242_metadata(meta);
            }
        let selected_face_set: BTreeSet<usize> = selection.selected_faces.iter().copied().collect();
        let mut selected_edge_set: BTreeSet<usize> = selection.selected_edges.iter().copied().collect();
        // If faces are selected, include their boundary edges in the 1D export so
        // wire entities are preserved alongside surface subsets.
        if !selected_face_set.is_empty() {
            let tbrep = self.tbrep();
            let mut face_index = 0usize;
            for ts in &tbrep.tshapes {
                let topods::TShape::Solid(sd) = &**ts else { continue };
                for shell_sr in &sd.shells {
                    if let topods::TShape::Shell(shd) = &*tbrep.tshapes[shell_sr.index] {
                        for face_sr in &shd.faces {
                            if selected_face_set.contains(&face_index) {
                                if let topods::TShape::Face(fd) = &*tbrep.tshapes[face_sr.index] {
                                    if let topods::TShape::Wire(w_outer) = &*tbrep.tshapes[fd.outer_wire.index] {
                                        selected_edge_set.extend(w_outer.edges.iter().map(|sr| sr.index));
                                    }
                                    for inner_sr in &fd.inner_wires {
                                        if let topods::TShape::Wire(w_inner) = &*tbrep.tshapes[inner_sr.index] {
                                            selected_edge_set.extend(w_inner.edges.iter().map(|sr| sr.index));
                                        }
                                    }
                                }
                            }
                            face_index += 1;
                        }
                    }
                }
            }
        }
        let export_all = selected_face_set.is_empty() && selected_edge_set.is_empty();

        // Check if this is a compound / compsolid by inspecting topods tshapes
        let tbrep = self.tbrep();
        let is_compound = tbrep.tshapes.iter().any(|ts| matches!(&**ts, topods::TShape::Compound(_)));
        let is_compsolid = tbrep.tshapes.iter().any(|ts| matches!(&**ts, topods::TShape::CompSolid(_)));

        let mut face_items = Vec::new();
        let mut solid_items = Vec::new();
        let mut compound_items = Vec::new();
        let mut compsolid_items = Vec::new();
        let mut shell_face_groups: Vec<Vec<u64>> = Vec::new();
        let mut has_triangle_fallback = false;
        // Map from face_index -> list of STEP ADVANCED_FACE ids (for color assignment)
        let mut face_step_ids: Vec<(usize, Vec<u64>)> = Vec::new();

        // Pre-collect solid/shell/face data from topods to avoid borrow conflicts.
        struct TopodsSolidData {
            solid_idx: usize,
            shells: Vec<TopodsShellData>,
        }
        struct TopodsShellData {
            sr_index: usize,
            face_indices: Vec<usize>,
        }
        let mut topods_solids: Vec<TopodsSolidData> = Vec::new();
        unsafe {
            let tbrep = &*self.topods_brep;
            let mut solid_idx = 0usize;
            for ts in tbrep.tshapes.iter() {
                let topods::TShape::Solid(sd) = &**ts else { continue };
                let mut shells = Vec::new();
                for shell_sr in &sd.shells {
                    if let topods::TShape::Shell(shd) = &*tbrep.tshapes[shell_sr.index] {
                        let face_indices: Vec<usize> = shd.faces.iter().map(|sr| sr.index).collect();
                        shells.push(TopodsShellData { sr_index: shell_sr.index, face_indices });
                    }
                }
                topods_solids.push(TopodsSolidData { solid_idx, shells });
                solid_idx += 1;
            }
        }

        let mut face_index = 0usize;
        for solid_data in &topods_solids {
            // First pass: collect faces for each shell.
            let mut shell_infos: Vec<(Vec<u64>, bool)> = Vec::new(); // (face_ids, is_closed)
            for shell_data in &solid_data.shells {
                let mut shell_faces = Vec::new();
                for &face_tshape_idx in &shell_data.face_indices {
                    if export_all || selected_face_set.contains(&face_index) {
                        let export =
                            self.write_face_topods(face_tshape_idx);
                        if export.used_triangle_fallback {
                            has_triangle_fallback = true;
                        }
                        face_step_ids.push((face_index, export.face_ids.clone()));
                        face_items.extend(export.face_ids.iter().copied());
                        shell_faces.extend(export.face_ids);
                    }
                    face_index += 1;
                }
                if !shell_faces.is_empty() {
                    let is_closed = unsafe {
                        let tbrep = &*self.topods_brep;
                        shell_is_closed_topods(tbrep, shell_data.sr_index)
                    };
                    shell_infos.push((shell_faces, is_closed));
                }
            }

            // Second pass: write topology 锟?MANIFOLD_SOLID_BREP for
            // single-shell solids, BREP_WITH_VOIDS for multi-shell solids
            // (outer shell + void shells), matching OCCT's export behaviour.
            let num_exported_solids = if export_all {
                let closed: Vec<&Vec<u64>> = shell_infos
                    .iter()
                    .filter(|(_, is_closed)| *is_closed)
                    .map(|(faces, _)| faces)
                    .collect();
                if closed.len() == 1 {
                    let shell_id = self.closed_shell(
                        &format!("closed_shell_{}_0", solid_data.solid_idx),
                        closed[0],
                    );
                    let solid_id = self.manifold_solid_brep(
                        &format!("solid_{}", solid_data.solid_idx),
                        shell_id,
                    );
                    solid_items.push(solid_id);
                } else if closed.len() > 1 {
                    let outer_id = self.closed_shell(
                        &format!("closed_shell_{}_0", solid_data.solid_idx),
                        closed[0],
                    );
                    let void_ids: Vec<u64> = closed[1..]
                        .iter()
                        .enumerate()
                        .map(|(i, faces)| {
                            self.closed_shell(
                                &format!("void_shell_{}_{}", solid_data.solid_idx, i),
                                faces,
                            )
                        })
                        .collect();
                    let solid_id = self.brep_with_voids(
                        &format!("solid_{}", solid_data.solid_idx),
                        outer_id,
                        &void_ids,
                    );
                    solid_items.push(solid_id);
                }
                closed.len()
            } else {
                0
            };

            // Shells that were NOT exported as MANIFOLD_SOLID_BREP / BREP_WITH_VOIDS
            // go to the SHELL_BASED_SURFACE_MODEL path (selected-face export, open shells).
            for (si, (faces, _is_closed)) in shell_infos.iter().enumerate() {
                if !faces.is_empty() && si >= num_exported_solids {
                    shell_face_groups.push(faces.clone());
                }
            }
        }

        // Handle compound / compsolid structure via topods inspection
        let has_compound = { let tbrep = self.tbrep(); tbrep.tshapes.iter().any(|ts| matches!(&**ts, topods::TShape::Compound(_))) };
        let has_compsolid = { let tbrep = self.tbrep(); tbrep.tshapes.iter().any(|ts| matches!(&**ts, topods::TShape::CompSolid(_))) };
        if has_compound && !solid_items.is_empty() {
            compound_items.push(self.compound("compound", &solid_items));
        }
        if has_compsolid && !solid_items.is_empty() && compound_items.is_empty() {
            compsolid_items.push(self.compsolid("compsolid", &solid_items));
        }

        if has_triangle_fallback {
            // Keep manifold solids for shells that were exported analytically.
            // A local triangle fallback should not force a full model downgrade.
        }

        // Face-boundary edges are already represented through face topology.
        // We keep this set to avoid duplicating them in wireframe export,
        // while still exporting standalone 1D curves.
        let mut face_edge_set: BTreeSet<usize> = BTreeSet::new();
        unsafe {
            let tbrep = &*self.topods_brep;
            for ts in &tbrep.tshapes {
                if let topods::TShape::Solid(sd) = &**ts {
                    for shell_sr in &sd.shells {
                        if let topods::TShape::Shell(shd) = &*tbrep.tshapes[shell_sr.index] {
                            for face_sr in &shd.faces {
                                if let topods::TShape::Face(fd) = &*tbrep.tshapes[face_sr.index] {
                                    if let topods::TShape::Wire(wd) = &*tbrep.tshapes[fd.outer_wire.index] {
                                        face_edge_set.extend(wd.edges.iter().map(|sr| sr.index));
                                    }
                                    for inner_sr in &fd.inner_wires {
                                        if let topods::TShape::Wire(wd) = &*tbrep.tshapes[inner_sr.index] {
                                            face_edge_set.extend(wd.edges.iter().map(|sr| sr.index));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Export standalone edge geometry (1D items not belonging to any face)
        // into wireframe representations. Use geometric curve entities instead
        // of topological EDGE_CURVE entities to align with OCCT-style
        // GEOMETRIC_CURVE_SET usage.
        let mut edge_items = Vec::new();
        let collect_standalone = export_all && (self.export_standalone_wire_overlay || !selected_edge_set.is_empty());
        unsafe {
            let tbrep = &*self.topods_brep;
            for (ti, ts) in tbrep.tshapes.iter().enumerate() {
                let topods::TShape::Edge(ed) = &**ts else { continue };
                let edge_idx = ti;
                if collect_standalone && face_edge_set.contains(&edge_idx) {
                    continue;
                }
                if collect_standalone {
                    if !ed.pcurves.is_empty() {
                        // Keep wireframe export focused on standalone 3D curves.
                        // SURFACE_CURVE + PCURVE-bound edges belong to face topology.
                        continue;
                    }
                    let start = tbrep.tshapes.get(ed.first.index).and_then(|ts| {
                        if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
                    });
                    let end = tbrep.tshapes.get(ed.last.index).and_then(|ts| {
                        if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
                    });
                    let (Some(start), Some(end)) = (start, end) else { continue };
                    if (end - start).length_squared() <= 1e-20 {
                        // Degenerate standalone curve segments can invalidate
                        // downstream wireframe representation parsing.
                        continue;
                    }
                    edge_items.push(self.write_standalone_wire_curve_by_index_topods(edge_idx));
                } else if selected_edge_set.contains(&edge_idx) {
                    edge_items.push(self.write_standalone_wire_curve_by_index_topods(edge_idx));
                }
            }
        }

        let mut standalone_point_items = Vec::new();
        if export_all && self.export_standalone_wire_overlay {
            let mut used_vertices: BTreeSet<usize> = BTreeSet::new();
            unsafe {
                let tbrep = &*self.topods_brep;
                for ts in &tbrep.tshapes {
                    if let topods::TShape::Edge(ed) = &**ts {
                        used_vertices.insert(ed.first.index);
                        used_vertices.insert(ed.last.index);
                    }
                }
                for (ti, ts) in tbrep.tshapes.iter().enumerate() {
                    if let topods::TShape::Vertex(v) = &**ts {
                        if !used_vertices.contains(&ti) {
                            standalone_point_items.push(
                                self.cartesian_point("", dvec3_to_array(v.point))
                            );
                        }
                    }
                }

                if standalone_point_items.is_empty() {
                    for ts in &tbrep.tshapes {
                        if let topods::TShape::Vertex(v) = &**ts {
                            standalone_point_items.push(
                                self.cartesian_point("", dvec3_to_array(v.point))
                            );
                        }
                        if standalone_point_items.len() >= 2 {
                            break;
                        }
                    }
                }
            }
        }

        // Application context strings depend on the selected protocol.
        let (ctx_name, proto_name, proto_year) = match self.protocol {
            StepProtocol::Ap214 => (
                "automotive_design",
                "automotive_design",
                2000_i32,
            ),
            StepProtocol::Ap242 => (
                "managed model based 3d engineering",
                "ap242_managed_model_based_3d_engineering_mim_lf",
                2014_i32,
            ),
        };
        let app_context = self.application_context(ctx_name);
        let _protocol = self.application_protocol_definition(
            "international standard",
            proto_name,
            proto_year,
            app_context,
        );
        let product_context = self.product_context("part definition", app_context, "mechanical");
        let product = self.product("rcad_export", "rcad_export", "", &[product_context]);
        let formation = self.product_definition_formation("", "", product);
        let definition_context =
            self.product_definition_context("part definition", app_context, "design");
        let definition = self.product_definition("", "", formation, definition_context);
        let product_shape = self.product_definition_shape("", "", definition);

        // B-rep vertex/surface coordinates are always kernel SI metres; `LENGTH_UNIT` uses si_metre.
        let length_unit = self.length_unit_meter();
        let angle_unit = self.plane_angle_unit_radian();
        let solid_angle_unit = self.solid_angle_unit_steradian();
        let write_tol = {
            let tbrep = self.tbrep();
            let max_vert = tbrep.tshapes.iter().filter_map(|ts| {
                if let topods::TShape::Vertex(vd) = &**ts { Some(vd.tolerance) } else { None }
            }).fold(0.0_f64, f64::max);
            let max_edge = tbrep.tshapes.iter().filter_map(|ts| {
                if let topods::TShape::Edge(ed) = &**ts { Some(ed.tolerance) } else { None }
            }).fold(0.0_f64, f64::max);
            let max_tol = max_vert.max(max_edge);
            if max_tol > 1e-12 { max_tol * 1.1 + 1e-7 } else { 1e-6 }
        };
        let uncertainty = self.uncertainty_measure_with_unit_value(length_unit, write_tol);
        let context = self.geometric_representation_context(
            3,
            uncertainty,
            &[length_unit, angle_unit, solid_angle_unit],
            "Context #1",
            "3D Context with UNIT and UNCERTAINTY",
        );

        let mut shell_model_items = Vec::new();
        if !face_items.is_empty() {
            for (i, shell_faces) in shell_face_groups.iter().enumerate() {
                if shell_faces.is_empty() {
                    continue;
                }
                let shell_id = self.open_shell(&format!("export_shell_{}", i), shell_faces);
                let model_id = self
                    .shell_based_surface_model(&format!("export_shell_model_{}", i), &[shell_id]);
                shell_model_items.push(model_id);
            }
        }

        if export_all {
            let origin = self.cartesian_point("asm_origin", [0.0, 0.0, 0.0]);
            let axis = self.direction("asm_axis", [0.0, 0.0, 1.0]);
            let ref_dir = self.direction("asm_ref", [1.0, 0.0, 0.0]);
            let root_axis = self.axis2_placement_3d("", origin, axis, ref_dir);

            let mut child_reps: Vec<(u64, u64)> = Vec::new();

            for (i, solid_id) in solid_items.iter().enumerate() {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, write_tol);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.advanced_brep_shape_representation(
                    &format!("strict_solid_rep_{}", i + 1),
                    &[child_axis, *solid_id],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            for (i, shell_model_id) in shell_model_items.iter().enumerate() {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, write_tol);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.manifold_surface_shape_representation(
                    &format!("strict_surface_rep_{}", i + 1),
                    &[child_axis, *shell_model_id],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            for (i, edge_item) in edge_items.iter().enumerate() {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let curve_set = self.geometric_curve_set("wireframe", std::slice::from_ref(edge_item));
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, write_tol);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.geometrically_bounded_wireframe_shape_representation(
                    &format!("strict_wire_rep_{}", i + 1),
                    &[child_axis, curve_set],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            if edge_items.len() > 1 {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let mut all_wire_items = edge_items.clone();
                all_wire_items.extend(standalone_point_items.iter().copied());
                let curve_set = self.geometric_curve_set("wireframe", &all_wire_items);
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, write_tol);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.geometrically_bounded_wireframe_shape_representation(
                    "strict_wire_rep_all",
                    &[child_axis, curve_set],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            let root_items = vec![root_axis];
            let root_rep = self.shape_representation("", &root_items, context);
            let _shape_def = self.shape_definition_representation(product_shape, root_rep);
            let _root_cat = self.product_related_product_category("part", None, &[product]);

            for (i, (child_rep, child_axis)) in child_reps.iter().enumerate() {
                let child_product_context = self.product_context("", app_context, "mechanical");
                let child_product = self.product(
                    &format!("rcad_export.{}", i + 1),
                    &format!("rcad_export.{}", i + 1),
                    "",
                    &[child_product_context],
                );
                let child_formation = self.product_definition_formation("", "", child_product);
                let child_def_context = self.product_definition_context("", app_context, "design");
                let child_definition =
                    self.product_definition("design", "", child_formation, child_def_context);
                let child_shape = self.product_definition_shape("", "", child_definition);
                let _child_shape_def = self.shape_definition_representation(child_shape, *child_rep);

                let placement_nauo = self.next_assembly_usage_occurrence(
                    &(i + 1).to_string(),
                    "",
                    "",
                    definition,
                    child_definition,
                );
                let placement_shape = self.product_definition_shape(
                    "Placement",
                    "Placement of an item",
                    placement_nauo,
                );
                let tr = self.item_defined_transformation("", "", root_axis, *child_axis);
                let rep_rel = self.representation_relationship_with_transformation(
                    "",
                    "",
                    *child_rep,
                    root_rep,
                    tr,
                );
                let _cdsr = self.context_dependent_shape_representation(rep_rel, placement_shape);
                let _child_cat = self.product_related_product_category("part", None, &[child_product]);
            }
        } else {
            let mut brep_rep = None;
            if export_all && !compound_items.is_empty() {
                brep_rep = Some(
                    self.advanced_brep_shape_representation("rcad_export", &compound_items, context),
                );
            } else if export_all && !compsolid_items.is_empty() {
                brep_rep = Some(
                    self.advanced_brep_shape_representation("rcad_export", &compsolid_items, context),
                );
            } else if export_all && !solid_items.is_empty() {
                brep_rep = Some(
                    self.advanced_brep_shape_representation("rcad_export", &solid_items, context),
                );
            }

            let mut surface_rep = None;
            if !shell_model_items.is_empty() {
                surface_rep = Some(self.manifold_surface_shape_representation(
                    "rcad_export",
                    &shell_model_items,
                    context,
                ));
            }

            let include_wire_overlay = brep_rep.is_none()
                || !selected_edge_set.is_empty()
                || !selected_face_set.is_empty()
                || (self.export_standalone_wire_overlay && export_all && !edge_items.is_empty());
            let mut wire_rep = None;
            if include_wire_overlay && !edge_items.is_empty() {
                let curve_sets: Vec<u64> = vec![self.geometric_curve_set("wireframe", &edge_items)];
                wire_rep = Some(self.geometrically_bounded_wireframe_shape_representation(
                    "rcad_export",
                    &curve_sets,
                    context,
                ));
                // Style standalone wire curves so they stay visible in viewers
                // that otherwise render default white-on-white curves.
                let wire_color = Color {
                    r: 0.560_784_34,
                    g: 0.686_274_5,
                    b: 0.560_784_34,
                };
                for &curve_id in &edge_items {
                    self.write_curve_color(curve_id, wire_color, 0.1);
                }
            }

            // When exporting standalone 1D entities with full solids, bridge through
            // a root SHAPE_REPRESENTATION like hfss/OCCT to improve viewer compatibility.
            let use_root_bridge = export_all && wire_rep.is_some() && brep_rep.is_some();
            if use_root_bridge {
                let root_origin = self.cartesian_point("asm_origin", [0.0, 0.0, 0.0]);
                let root_axis_dir = self.direction("asm_axis", [0.0, 0.0, 1.0]);
                let root_ref_dir = self.direction("asm_ref", [1.0, 0.0, 0.0]);
                let root_axis = self.axis2_placement_3d("", root_origin, root_axis_dir, root_ref_dir);
                let root_rep = self.shape_representation("rcad_export", &[root_axis], context);
                let _shape_def = self.shape_definition_representation(product_shape, root_rep);

                if let Some(rep) = brep_rep {
                    self.shape_representation_relationship("", "", root_rep, rep);
                }
                if let Some(rep) = wire_rep {
                    self.shape_representation_relationship("", "", root_rep, rep);
                }
                if let Some(rep) = surface_rep {
                    self.shape_representation_relationship("", "", root_rep, rep);
                }
            } else {
                let primary_rep = brep_rep.or(surface_rep).or(wire_rep);
                let mut secondary_reps = Vec::new();
                if let Some(rep) = brep_rep && Some(rep) != primary_rep {
                    secondary_reps.push(rep);
                }
                if let Some(rep) = wire_rep && Some(rep) != primary_rep {
                    secondary_reps.push(rep);
                }
                if let Some(rep) = surface_rep && Some(rep) != primary_rep {
                    secondary_reps.push(rep);
                }

                let primary_rep =
                    primary_rep.unwrap_or_else(|| self.shape_representation("rcad_export", &[], context));
                let _shape_def = self.shape_definition_representation(product_shape, primary_rep);
                for rep in secondary_reps {
                    self.shape_representation_relationship("", "", primary_rep, rep);
                }
            }
        }

        // -- Color / presentation styling --
        if let Some(step_colors) = colors {
            for (fi, step_ids) in &face_step_ids {
                if let Some(color) = step_colors.color_for_face(*fi) {
                    for &face_id in step_ids {
                        self.write_face_color(face_id, color);
                    }
                }
            }
        } else if !face_step_ids.is_empty() {
            // Keep exports visible in viewers that default to white background
            // and don't auto-assign contrasting face colors.
            let fallback = Color {
                r: 0.560_784_34,
                g: 0.686_274_5,
                b: 0.560_784_34,
            };
            for (_, step_ids) in &face_step_ids {
                for &face_id in step_ids {
                    self.write_face_color(face_id, fallback);
                }
            }
        }

    }

    /// Write a compound structure to STEP, returning the list of STEP entity IDs.
    fn write_compound_structure(
        &mut self,
        compound: &flat::Compound,
        existing_solids: &[u64],
    ) -> Vec<u64> {
        let mut element_ids = Vec::new();

        // Use existing solids if available
        let mut solid_iter = existing_solids.iter();

        // Add solids
        for (i, (_label, _solid)) in compound.solids.iter().enumerate() {
            if let Some(&solid_id) = solid_iter.next() {
                element_ids.push(solid_id);
            } else {
                // Create a placeholder solid ID
                let id = self.manifold_solid_brep(&format!("compound_solid_{}", i), 0);
                element_ids.push(id);
            }
        }

        // Add comp_solids
        for (i, (_label, compsolid)) in compound.comp_solids.iter().enumerate() {
            let compsolid_id = self.write_compsolid_structure(compsolid, existing_solids);
            if compsolid_id.len() == 1 {
                element_ids.push(compsolid_id[0]);
            } else {
                let id = self.compsolid(&format!("compound_compsolid_{}", i), &compsolid_id);
                element_ids.push(id);
            }
        }

        // Add nested compounds recursively
        for (i, (_label, nested)) in compound.compounds.iter().enumerate() {
            let nested_items = self.write_compound_structure(nested, existing_solids);
            let compound_id = self.compound(&format!("nested_compound_{}", i), &nested_items);
            element_ids.push(compound_id);
        }

        // If we have multiple elements, wrap in a COMPOUND entity
        if element_ids.len() > 1 {
            let compound_id = self.compound("compound", &element_ids);
            vec![compound_id]
        } else {
            element_ids
        }
    }

    /// Write a compsolid structure to STEP, returning the list of STEP entity IDs.
    fn write_compsolid_structure(
        &mut self,
        compsolid: &flat::CompSolid,
        existing_solids: &[u64],
    ) -> Vec<u64> {
        let mut solid_ids = Vec::new();
        let mut solid_iter = existing_solids.iter();

        for (i, _solid) in compsolid.solids.iter().enumerate() {
            if let Some(&solid_id) = solid_iter.next() {
                solid_ids.push(solid_id);
            } else {
                let id = self.manifold_solid_brep(&format!("compsolid_solid_{}", i), 0);
                solid_ids.push(id);
            }
        }

        if solid_ids.len() > 1 {
            let compsolid_id = self.compsolid("compsolid", &solid_ids);
            vec![compsolid_id]
        } else {
            solid_ids
        }
    }

    /// Topods-native face writer: reads face data from tshapes instead of FlatBRep.
    fn write_face_topods(&mut self, face_tshape_idx: usize) -> FaceExportResult {
        // Phase 1: extract all data from tshapes (no self calls).
        let (face_surface, oriented_edges, loop_points, origin_point, normal, face_ref_arr, seam_edge_indices, inner_wires_data)
            = {
            let tbrep = self.tbrep();
            let topods::TShape::Face(fd) = &*tbrep.tshapes[face_tshape_idx] else {
                return FaceExportResult { face_ids: vec![], used_triangle_fallback: false };
            };
            let face_surface = fd.surface.clone();

            // Degenerate face — no surface and few edges
            if is_degenerate_face_wire_topods(tbrep, face_tshape_idx) {
                return FaceExportResult {
                    face_ids: vec![],
                    used_triangle_fallback: false,
                };
            }

            let oriented_edges = oriented_face_edges_topods(tbrep, face_tshape_idx);
            if oriented_edges.is_empty() && face_surface.is_some() {
                return FaceExportResult {
                    face_ids: vec![],
                    used_triangle_fallback: false,
                };
            }

            // Loop points from vertex tshape points
            let loop_points: Vec<glam::DVec3> = oriented_edges
                .iter()
                .filter_map(|edge| {
                    tbrep.tshapes.get(edge.start).and_then(|ts| {
                        if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
                    })
                })
                .collect();
            let origin_point = loop_points.first().copied().unwrap_or(glam::DVec3::ZERO);

            let normal = compute_face_normal(&loop_points)
                .or_else(|| surface_normal(face_surface.clone()))
                .map(dvec3_to_array)
                .unwrap_or([0.0, 0.0, 1.0]);
            let face_ref_arr = orthogonal_dir(normal);

            // Seam edges by edge tshape index
            let seam_edge_indices = detect_seam_edge_indices_topods(tbrep, face_tshape_idx);

            // Extract inner wire data (Vec of Vec<OrientedEdgeExport>)
            let inner_wires_data: Vec<Vec<OrientedEdgeExport>> = fd.inner_wires.iter().filter_map(|iw_sr| {
                let topods::TShape::Wire(iwd) = &*tbrep.tshapes[iw_sr.index] else { return None };
                let oriented: Vec<OrientedEdgeExport> = iwd.edges.iter().filter_map(|sr| {
                    let topods::TShape::Edge(ed) = &*tbrep.tshapes[sr.index] else { return None };
                    let (start, end) = if sr.orientation.is_forward() {
                        (ed.first.index, ed.last.index)
                    } else {
                        (ed.last.index, ed.first.index)
                    };
                    Some(OrientedEdgeExport { edge_idx: sr.index, start, end, forward: sr.orientation.is_forward() })
                }).collect();
                if oriented.is_empty() { None } else { Some(oriented) }
            }).collect();

            (face_surface, oriented_edges, loop_points, origin_point, normal, face_ref_arr, seam_edge_indices, inner_wires_data)
        };

        // Phase 2: write operations (self method calls, no tbrep borrow active).
        let surface = if face_surface.is_some() {
            self.get_or_write_surface_id_topods(face_tshape_idx)
        } else {
            let fallback_placement = {
                let origin = self.cartesian_point("face_origin", dvec3_to_array(origin_point));
                let axis = self.direction("face_normal", normal);
                let ref_dir = self.direction("face_ref", face_ref_arr);
                Some(self.axis2_placement_3d("face_axis", origin, axis, ref_dir))
            };
            self.write_surface(None, fallback_placement)
        };

        let face_orientation = {
            let tbrep = self.tbrep();
            face_orientation_for_surface_topods(tbrep, &loop_points, face_tshape_idx)
        };

        // Separate degenerate edges (self-loop) from normal edges
        let mut degen_verts: Vec<usize> = Vec::new();
        let mut edge_entries: Vec<(usize, usize, usize, u64, bool)> = Vec::new();
        for edge in &oriented_edges {
            if edge.start == edge.end {
                degen_verts.push(edge.start);
                continue;
            }
            let edge_curve = if seam_edge_indices.contains(&edge.edge_idx) {
                self.write_seam_edge_curve_topods(edge.edge_idx, face_surface.clone())
            } else {
                // write_edge_curve_by_index_topods handles pcurve seeding internally
                self.write_edge_curve_by_index_topods(edge.edge_idx)
            };
            edge_entries.push((edge.edge_idx, edge.start, edge.end, edge_curve, edge.forward));
        }

        let mut oriented_entries: Vec<(u64, bool)> = edge_entries
            .iter()
            .map(|(_, _, _, curve, forward)| (*curve, *forward))
            .collect();

        // OCCT-style cyclic ordering on periodic cylinder/cone side faces
        if matches!(face_surface.as_ref(), Some(Surface3::Cylinder(_)) | Some(Surface3::Cone(_)))
            && seam_edge_indices.len() == 1
            && edge_entries.len() == 4
        {
            let seam_idx = *seam_edge_indices.iter().next().unwrap_or(&usize::MAX);
            if seam_idx != usize::MAX {
                let seam_curve = edge_entries
                    .iter()
                    .find_map(|(idx, _, _, curve, _)| if *idx == seam_idx { Some(*curve) } else { None });
                let circle_entries: Vec<(usize, usize, usize, u64, bool)> = edge_entries
                    .iter()
                    .copied()
                    .filter(|(idx, _, _, _, _)| *idx != seam_idx)
                    .collect();
                if let (Some(seam_curve), 2) = (seam_curve, circle_entries.len()) {
                    let axis = face_surface
                        .as_ref()
                        .and_then(|surf| match surf {
                            Surface3::Cylinder(cyl) => Some(canonicalize_axis_sign(cyl.axis.normalize_or_zero())),
                            Surface3::Cone(cone) => Some(canonicalize_axis_sign(cone.axis_dir())),
                            _ => None,
                        })
                        .filter(|a: &glam::DVec3| a.length_squared() > 1e-18)
                        .unwrap_or_else(|| canonicalize_axis_sign(glam::DVec3::Z));
                    let z_of = |tbrep: &topods::BRep, vid: usize| -> f64 {
                        tbrep.tshapes.get(vid).and_then(|ts| {
                            if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
                        }).map(|p| p.dot(axis)).unwrap_or(0.0)
                    };
                    let tbrep_local = self.tbrep();
                    let (c0, c1) = (circle_entries[0], circle_entries[1]);
                    let top_circle = if z_of(tbrep_local, c0.1) >= z_of(tbrep_local, c1.1) { c0 } else { c1 };
                    let bottom_circle = if z_of(tbrep_local, c0.1) < z_of(tbrep_local, c1.1) { c0 } else { c1 };
                    oriented_entries = vec![
                        (top_circle.3, false),
                        (seam_curve, false),
                        (bottom_circle.3, true),
                        (seam_curve, true),
                    ];
                }
            }
        }

        let mut oriented_ids = Vec::with_capacity(oriented_entries.len());
        for (curve, orientation) in oriented_entries {
            oriented_ids.push(self.oriented_edge("face_edge", curve, orientation));
        }

        let edge_loop = self.edge_loop("outer_loop", &oriented_ids);
        let mut bounds = vec![self.face_bound("outer_bound", edge_loop, true)];

        for &dvi in &degen_verts {
            let vp = self.vertex_point_by_tshape_idx(dvi);
            let vl = self.vertex_loop("degen_loop", vp);
            bounds.push(self.face_bound("degen_bound", vl, true));
        }

        // SplitCommonVertex: collect outer wire vertices
        let outer_verts: HashSet<usize> = oriented_edges
            .iter()
            .flat_map(|e| [e.start, e.end])
            .collect();
        let mut all_wire_verts = outer_verts;

        // Inner wires (data already extracted in Phase 1)
        for (ii, inner_oriented) in inner_wires_data.iter().enumerate() {
            // Clear vertex_point_ids for shared vertices
            for oe in inner_oriented {
                if !all_wire_verts.insert(oe.start) {
                    self.vertex_point_ids.remove(&oe.start);
                }
                if !all_wire_verts.insert(oe.end) {
                    self.vertex_point_ids.remove(&oe.end);
                }
            }

            let mut inner_entries: Vec<(u64, bool)> = Vec::with_capacity(inner_oriented.len());
            for edge in inner_oriented {
                let curve = self.write_edge_curve_by_index_topods(edge.edge_idx);
                inner_entries.push((curve, edge.forward));
            }
            let mut inner_ids = Vec::with_capacity(inner_entries.len());
            for (curve, orientation) in inner_entries {
                inner_ids.push(self.oriented_edge("face_edge", curve, orientation));
            }
            let inner_loop = self.edge_loop(&format!("inner_loop_{ii}"), &inner_ids);
            bounds.push(self.face_bound(&format!("inner_bound_{ii}"), inner_loop, true));
        }

        FaceExportResult {
            face_ids: vec![self.advanced_face("face", &bounds, surface, face_orientation)],
            used_triangle_fallback: false,
        }
    }

    fn write_face(
        &mut self,
        brep: &BRep,
        face: &Face,
        face_surface: Option<Surface3>,
        face_surface_idx: Option<usize>,
    ) -> FaceExportResult {
        // Only fall back to triangles when there is no analytic surface AND the
        // wire is degenerate (cannot form a valid edge loop).
        if face_surface.is_none() && is_degenerate_face_wire(brep, face) {
            return FaceExportResult {
                face_ids: self.write_triangle_faces(brep, face),
                used_triangle_fallback: true,
            };
        }

        let oriented_edges = oriented_face_edges(brep, face);
        if oriented_edges.is_empty() && face_surface.is_some() {
            // Seam face with no usable edge loop ??fall back to triangles.
            return FaceExportResult {
                face_ids: self.write_triangle_faces(brep, face),
                used_triangle_fallback: true,
            };
        }

        let loop_points: Vec<glam::DVec3> = oriented_edges
            .iter()
            .filter_map(|edge| brep.vertices.get(edge.start).map(|v| v.point))
            .collect();
        let origin_point = loop_points.first().copied().unwrap_or(glam::DVec3::ZERO);

        // For seam faces on closed surfaces, the loop points may be collinear,
        // so compute_face_normal can fail. Use the surface's own axis instead.
        let normal = compute_face_normal(&loop_points)
            .or_else(|| surface_normal(face_surface.clone()))
            .map(dvec3_to_array)
            .unwrap_or([0.0, 0.0, 1.0]);

        let face_ref_arr = orthogonal_dir(normal);
        let surface = if let Some(surface_idx) = face_surface_idx {
            self.get_or_write_surface_id(brep, surface_idx)
        } else {
            let fallback_placement = if face_surface.is_none() {
                let origin = self.cartesian_point("face_origin", dvec3_to_array(origin_point));
                let axis = self.direction("face_normal", normal);
                let ref_dir = self.direction("face_ref", face_ref_arr);
                Some(self.axis2_placement_3d("face_axis", origin, axis, ref_dir))
            } else {
                None
            };
            self.write_surface(face_surface.clone(), fallback_placement)
        };
        let face_orientation = face_orientation_for_surface(
            brep,
            &loop_points,
            face,
            face_surface.as_ref(),
        );

        // Detect seam edges: same edge_idx appearing multiple times
        let seam_edge_indices = detect_seam_edge_indices(face);

        // Separate degenerate edges (self-loop, start==end) from normal edges.
        // OCCT writes degenerate edges as FACE_BOUND + VERTEX_LOOP instead of
        // EDGE_CURVE, so skip them from the edge loop to avoid breaking shell closure.
        let mut degen_verts: Vec<usize> = Vec::new();
        let mut edge_entries: Vec<(usize, usize, usize, u64, bool)> = Vec::new();
        for edge in &oriented_edges {
            if edge.start == edge.end {
                degen_verts.push(edge.start);
                continue;
            }
            let edge_curve = if seam_edge_indices.contains(&edge.edge_idx) {
                // Write a proper SEAM_CURVE (or SURFACE_CURVE with 2+ pcurves)
                // for any edge that appears multiple times in the same face wire.
                // The previous logic branched on has_surface_parametrics and
                // wrote SURFACE_CURVE when pcurves existed 锟?this prevented
                // SEAM_CURVE output for synthetic seams (e.g. sphere-box union
                // analytic builder) that explicitly provide sphere pcurves.
                self.write_seam_edge_curve_cached(brep, edge.edge_idx, face_surface.clone())
            } else {
                self.seed_edge_surface_curve_from_face_frame(
                    brep,
                    edge.edge_idx,
                    face_surface_idx,
                    surface,
                    origin_point,
                    normal,
                    face_ref_arr,
                );
                self.write_edge_curve_by_index_topods(edge.edge_idx)
            };
            edge_entries.push((edge.edge_idx, edge.start, edge.end, edge_curve, edge.forward));
        }

        let mut oriented_entries: Vec<(u64, bool)> = edge_entries
            .iter()
            .map(|(_, _, _, curve, forward)| (*curve, *forward))
            .collect();

        // OCCT-style cyclic ordering on periodic cylinder/cone side faces:
        // top circle (F) -> seam (F) -> bottom circle (T) -> seam (T)
        if matches!(face_surface.as_ref(), Some(Surface3::Cylinder(_)) | Some(Surface3::Cone(_)))
            && seam_edge_indices.len() == 1
            && edge_entries.len() == 4
        {
            let seam_idx = *seam_edge_indices.iter().next().unwrap_or(&usize::MAX);
            if seam_idx != usize::MAX {
                let seam_curve = edge_entries
                    .iter()
                    .find_map(|(idx, _, _, curve, _)| if *idx == seam_idx { Some(*curve) } else { None });
                let circle_entries: Vec<(usize, usize, usize, u64, bool)> = edge_entries
                    .iter()
                    .copied()
                    .filter(|(idx, _, _, _, _)| *idx != seam_idx)
                    .collect();

                if let (Some(seam_curve), 2) = (seam_curve, circle_entries.len()) {
                    let axis = face_surface
                        .as_ref()
                        .and_then(|surf| match surf {
                            Surface3::Cylinder(cyl) => Some(canonicalize_axis_sign(cyl.axis.normalize_or_zero())),
                            Surface3::Cone(cone) => Some(canonicalize_axis_sign(cone.axis_dir())),
                            _ => None,
                        })
                        .filter(|axis: &glam::DVec3| axis.length_squared() > 1e-18)
                        .unwrap_or_else(|| canonicalize_axis_sign(glam::DVec3::Z));
                    let z_of = |vid: usize| -> f64 {
                        brep.vertices
                            .get(vid)
                            .map(|v| v.point.dot(axis))
                            .unwrap_or(0.0)
                    };
                    let (c0, c1) = (circle_entries[0], circle_entries[1]);
                    let top_circle = if z_of(c0.1) >= z_of(c1.1) { c0 } else { c1 };
                    let bottom_circle = if z_of(c0.1) < z_of(c1.1) { c0 } else { c1 };

                    oriented_entries = vec![
                        (top_circle.3, false),
                        (seam_curve, false),
                        (bottom_circle.3, true),
                        (seam_curve, true),
                    ];
                }
            }
        }

        let mut oriented_ids = Vec::with_capacity(oriented_entries.len());
        for (curve, orientation) in oriented_entries {
            oriented_ids.push(self.oriented_edge("face_edge", curve, orientation));
        }

        let edge_loop = self.edge_loop("outer_loop", &oriented_ids);
        let mut bounds = vec![self.face_bound("outer_bound", edge_loop, true)];

        // OCCT-aligned: degenerate edges (self-loops) are written as FACE_BOUND
        // with VERTEX_LOOP, not as EDGE_CURVE, to keep shell closure intact.
        for &dvi in &degen_verts {
            let vp = self.vertex_point_by_index(brep, dvi);
            let vl = self.vertex_loop("degen_loop", vp);
            bounds.push(self.face_bound("degen_bound", vl, true));
        }

        // 鈹€鈹€ SplitCommonVertex 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        // STEP requires vertices belonging to different FACE_BOUND entities
        // (outer wire vs inner wires) to be distinct VERTEX_POINT entities.
        // Collect outer wire vertices and clear them from the cache before
        // processing inner wires so that shared geometric vertices produce
        // separate VERTEX_POINT entities.  Analogous to OCCT's
        // `ShapeFix_SplitCommonVertex`.
        let outer_verts: HashSet<usize> = oriented_edges
            .iter()
            .flat_map(|e| [e.start, e.end])
            .collect();
        let mut all_wire_verts = outer_verts;

        // Inner wires (holes): for each inner wire, create an edge loop and face_bound.
        for (ii, inner_wire) in face.inner_wires.iter().enumerate() {
            let inner_oriented: Vec<OrientedEdgeExport> = inner_wire.edges.iter().filter_map(|we| {
                let edge = brep.edges.get(we.idx)?;
                let (start, end) = if we.forward { (edge.start, edge.end) } else { (edge.end, edge.start) };
                Some(OrientedEdgeExport { edge_idx: we.idx, start, end, forward: we.forward })
            }).collect();
            if inner_oriented.is_empty() { continue; }

            // Clear vertex_point_ids for vertices shared with any previously
            // written wire of this face (outer or earlier inner wires) so that
            // write_edge_curve_by_index creates unique VERTEX_POINT entities.
            for oe in &inner_oriented {
                if !all_wire_verts.insert(oe.start) {
                    self.vertex_point_ids.remove(&oe.start);
                }
                if !all_wire_verts.insert(oe.end) {
                    self.vertex_point_ids.remove(&oe.end);
                }
            }

            let mut inner_entries: Vec<(u64, bool)> = Vec::with_capacity(inner_oriented.len());
            for edge in &inner_oriented {
                self.seed_edge_surface_curve_from_face_frame(
                    brep,
                    edge.edge_idx,
                    face_surface_idx,
                    surface,
                    origin_point,
                    normal,
                    face_ref_arr,
                );
                let curve = self.write_edge_curve_by_index_topods(edge.edge_idx);
                let orientation = edge.forward;
                inner_entries.push((curve, orientation));
            }
            let mut inner_ids = Vec::with_capacity(inner_entries.len());
            for (curve, orientation) in inner_entries {
                inner_ids.push(self.oriented_edge("face_edge", curve, orientation));
            }
            let inner_loop = self.edge_loop(&format!("inner_loop_{ii}"), &inner_ids);
            bounds.push(self.face_bound(&format!("inner_bound_{ii}"), inner_loop, true));
        }
        FaceExportResult {
            face_ids: vec![self.advanced_face(
                "face",
                &bounds,
                surface,
                face_orientation,
            )],
            used_triangle_fallback: false,
        }
    }

    fn write_seam_edge_curve_cached(
        &mut self,
        _brep: &BRep,
        edge_idx: usize,
        face_surface: Option<Surface3>,
    ) -> u64 {
        if let Some(existing) = self.seam_edge_curve_ids.get(&edge_idx) {
            return *existing;
        }
        let edge_id = self.write_seam_edge_curve_topods(edge_idx, face_surface);
        self.seam_edge_curve_ids.insert(edge_idx, edge_id);
        edge_id
    }

    fn write_triangle_faces(&mut self, brep: &BRep, face: &Face) -> Vec<u64> {
        let mut faces = Vec::new();
        for tri in &face.triangles {
            let Some(a) = brep.vertices.get(tri[0]).map(|v| v.point) else {
                continue;
            };
            let Some(b) = brep.vertices.get(tri[1]).map(|v| v.point) else {
                continue;
            };
            let Some(c) = brep.vertices.get(tri[2]).map(|v| v.point) else {
                continue;
            };

            let n = (b - a).cross(c - a).normalize_or_zero();
            if n.length_squared() < 1e-12 {
                continue;
            }

            let origin = self.cartesian_point("tri_origin", dvec3_to_array(a));
            let axis = self.direction("tri_normal", dvec3_to_array(n));
            let ref_dir = self.direction("tri_ref", orthogonal_dir(dvec3_to_array(n)));
            let placement = self.axis2_placement_3d("tri_axis", origin, axis, ref_dir);
            let plane = self.plane("tri_plane", placement);

            let e0 = self.write_edge_curve_from_points(a, b);
            let e1 = self.write_edge_curve_from_points(b, c);
            let e2 = self.write_edge_curve_from_points(c, a);
            let o0 = self.oriented_edge("tri_edge", e0, true);
            let o1 = self.oriented_edge("tri_edge", e1, true);
            let o2 = self.oriented_edge("tri_edge", e2, true);
            let loop_id = self.edge_loop("tri_loop", &[o0, o1, o2]);
            let bound = self.face_outer_bound("tri_bound", loop_id, true);
            faces.push(self.advanced_face("tri_face", &[bound], plane, true));
        }
        faces
    }

    fn write_edge_curve_from_points(&mut self, a: glam::DVec3, b: glam::DVec3) -> u64 {
        let p0 = self.cartesian_point("tri_p0", dvec3_to_array(a));
        let p1 = self.cartesian_point("tri_p1", dvec3_to_array(b));
        let v0 = self.vertex_point("tri_v0", p0);
        let v1 = self.vertex_point("tri_v1", p1);
        let delta = dvec3_to_array(b - a);
        let dir = self.direction("tri_dir", normalize(delta));
        let vec = self.vector("tri_vec", dir, vector_length(delta).max(1e-9));
        let line = self.line("tri_line", p0, vec);
        self.edge_curve("tri_edge", v0, v1, line, true)
    }

    fn seed_edge_surface_curve_from_face_frame(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        face_surface_idx: Option<usize>,
        surface_id: u64,
        face_origin: glam::DVec3,
        face_normal: [f64; 3],
        face_ref: [f64; 3],
    ) {
        if self.edge_geometry_ids.contains_key(&edge_idx) {
            return;
        }
        if brep
            .geom
            .edge_pcurves
            .get(edge_idx)
            .is_some_and(|pcs| !pcs.is_empty())
        {
            return;
        }

        let Some(edge) = brep.edges.get(edge_idx) else {
            return;
        };
        let start_point = brep
            .vertices
            .get(edge.start)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let end_point = brep
            .vertices
            .get(edge.end)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);

        let u_axis = glam::DVec3::from_array(normalize(face_ref));
        let n_axis = glam::DVec3::from_array(normalize(face_normal));
        let mut v_axis = n_axis.cross(u_axis);
        if v_axis.length_squared() < 1e-18 {
            v_axis = any_perpendicular_dvec3(n_axis);
        } else {
            v_axis = v_axis.normalize_or_zero();
        }

        let current_surface = face_surface_idx.and_then(|idx| brep.geom.surfaces.get(idx).cloned());
        let curve2d = match current_surface.as_ref() {
            Some(Surface3::Plane(plane)) => synthesize_plane_pcurve_for_edge(brep, edge_idx, plane),
            Some(Surface3::Cylinder(cyl)) => synthesize_cylinder_pcurve_for_edge(brep, edge_idx, cyl),
            _ => synthesize_edge_curve2d_on_face_frame(
                brep,
                edge_idx,
                face_origin,
                u_axis,
                n_axis,
            )
            .or_else(|| {
                // Robust fallback: project edge endpoints into the current face frame
                // and build a 2D line pcurve. This keeps strict export aligned with
                // OCCT-style expectation that face edges carry parametric curves.
                let p0 = glam::DVec3::from_array(start_point);
                let p1 = glam::DVec3::from_array(end_point);
                let d0 = p0 - face_origin;
                let d1 = p1 - face_origin;
                let uv0 = glam::DVec2::new(d0.dot(u_axis), d0.dot(v_axis));
                let uv1 = glam::DVec2::new(d1.dot(u_axis), d1.dot(v_axis));
                let dir = uv1 - uv0;
                if dir.length_squared() < 1e-24 {
                    None
                } else {
                    Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                        origin: uv0,
                        direction: dir,
                    }))
                }
            }),
        };
        let Some(mut curve2d) = curve2d else {
            return;
        };
        let curve2d_was_line = matches!(curve2d, Curve2d::Line(_));
        let curve2d_was_circle = matches!(curve2d, Curve2d::Circle(_));
        let should_promote_current_plane_line = curve2d_was_line
            && should_promote_plane_line_pcurve(brep, edge_idx);
        let should_promote_current_cylinder_line = matches!(current_surface.as_ref(), Some(Surface3::Cylinder(cyl))
            if should_promote_cylinder_line_pcurve(cyl, glam::DVec3::from_array(start_point)));

        if matches!(current_surface.as_ref(), Some(Surface3::Plane(_)))
            && should_promote_current_plane_line
        {
            curve2d = plane_line_pcurve_as_bspline(&curve2d);
        }

        if should_promote_current_cylinder_line
            && curve2d_was_line
            && count_cylinder_face_occurrences_for_edge(brep, edge_idx) >= 2
        {
            curve2d = cylinder_line_pcurve_as_bspline(&curve2d);
        }

        let basis_curve_3d = self.write_basis_curve_for_edge(brep, edge_idx, start_point, end_point);
        let param_curve_id = self.write_curve2d(Some(curve2d.clone()));
        let def_rep = self.definitional_representation(param_curve_id);
        let pcurve_id = self.pcurve(surface_id, def_rep);
        let mut pcurve_ids = vec![pcurve_id];

        if curve2d_was_line
            && count_plane_face_occurrences_for_line_edge(brep, edge_idx) >= 2
        {
            let base_plane = rcad_kernel::geom::Plane {
                origin: face_origin,
                normal: n_axis,
            };
            let extra_surface = find_peer_plane_surface_for_line_edge(brep, edge_idx, base_plane)
                .or(Some((None, base_plane)));

            if let Some((peer_surface_idx, peer_plane)) = extra_surface
                && let Some(extra_curve2d) =
                    synthesize_plane_pcurve_for_edge(brep, edge_idx, &peer_plane)
            {
                let extra_curve2d = if should_promote_current_plane_line {
                    plane_line_pcurve_as_bspline(&extra_curve2d)
                } else {
                    extra_curve2d
                };
                let peer_surface_id = if let Some(idx) = peer_surface_idx {
                    self.get_or_write_surface_id(brep, idx)
                } else {
                    // Reuse current face surface id to avoid emitting synthetic
                    // duplicate PLANE entities for strict dual-PCURVE fallback.
                    surface_id
                };
                let extra_param = self.write_curve2d(Some(extra_curve2d));
                let extra_def = self.definitional_representation(extra_param);
                let extra_pc = self.pcurve(peer_surface_id, extra_def);
                pcurve_ids.push(extra_pc);
            }
        }

        if matches!(current_surface.as_ref(), Some(Surface3::Plane(_)))
            && curve2d_was_circle
            && let Some((cyl_surface_idx, cyl)) =
                find_cylinder_surface_for_edge_excluding(brep, edge_idx, face_surface_idx)
            && let Some(cyl_curve2d) = synthesize_cylinder_pcurve_for_edge(brep, edge_idx, &cyl)
        {
            let cyl_surface_id = self.get_or_write_surface_id(brep, cyl_surface_idx);
            let cyl_param = self.write_curve2d(Some(cyl_curve2d));
            let cyl_def = self.definitional_representation(cyl_param);
            let cyl_pc = self.pcurve(cyl_surface_id, cyl_def);
            pcurve_ids.push(cyl_pc);
        }

        if matches!(current_surface.as_ref(), Some(Surface3::Plane(_)))
            && curve2d_was_line
            && let Some((cyl_surface_idx, cyl)) =
                find_cylinder_surface_for_edge_excluding(brep, edge_idx, face_surface_idx)
            && let Some(cyl_curve2d) = synthesize_cylinder_pcurve_for_edge(brep, edge_idx, &cyl)
        {
            let cyl_surface_id = self.get_or_write_surface_id(brep, cyl_surface_idx);
            let cyl_curve2d = if should_promote_cylinder_line_pcurve(&cyl, glam::DVec3::from_array(start_point)) {
                cylinder_line_pcurve_as_bspline(&cyl_curve2d)
            } else {
                cyl_curve2d
            };
            let cyl_param = self.write_curve2d(Some(cyl_curve2d));
            let cyl_def = self.definitional_representation(cyl_param);
            let cyl_pc = self.pcurve(cyl_surface_id, cyl_def);
            pcurve_ids.push(cyl_pc);
        }

        if matches!(current_surface.as_ref(), Some(Surface3::Cylinder(_)))
            && curve2d_was_line
            && let Some((peer_surface_idx, peer_cyl)) =
                find_peer_cylinder_surface_for_edge(brep, edge_idx, face_surface_idx)
            && let Some(peer_curve2d) = synthesize_cylinder_pcurve_for_edge(brep, edge_idx, &peer_cyl)
        {
            let peer_surface_id = self.get_or_write_surface_id(brep, peer_surface_idx);
            let peer_curve2d = if count_cylinder_face_occurrences_for_edge(brep, edge_idx) >= 2
                && should_promote_cylinder_line_pcurve(&peer_cyl, glam::DVec3::from_array(start_point)) {
                cylinder_line_pcurve_as_bspline(&peer_curve2d)
            } else {
                peer_curve2d
            };
            let peer_param = self.write_curve2d(Some(peer_curve2d));
            let peer_def = self.definitional_representation(peer_param);
            let peer_pc = self.pcurve(peer_surface_id, peer_def);
            pcurve_ids.push(peer_pc);
        }

        if matches!(current_surface.as_ref(), Some(Surface3::Cylinder(_)))
            && curve2d_was_line
            && let Some((plane_surface_idx, plane)) = find_plane_surface_for_edge(brep, edge_idx)
            && let Some(plane_curve2d) = synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane)
        {
            let plane_surface_id = self.get_or_write_surface_id(brep, plane_surface_idx);
            let plane_curve2d = if should_promote_plane_line_pcurve(brep, edge_idx) {
                plane_line_pcurve_as_bspline(&plane_curve2d)
            } else {
                plane_curve2d
            };
            let plane_param = self.write_curve2d(Some(plane_curve2d));
            let plane_def = self.definitional_representation(plane_param);
            let plane_pc = self.pcurve(plane_surface_id, plane_def);
            pcurve_ids.push(plane_pc);
        }

        let surface_curve = self.surface_curve(basis_curve_3d, &pcurve_ids);
        self.edge_geometry_ids.insert(edge_idx, surface_curve);
    }

    /// Write an EDGE_CURVE for a seam edge from topods data.
    /// Mirrors write_seam_edge_curve but reads edge and vertex data
    /// from self.tbrep() (topods) instead of a flat BRep.
    fn write_seam_edge_curve_topods(
        &mut self,
        edge_idx: usize,
        face_surface: Option<Surface3>,
    ) -> u64 {
        let ed = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(edge_idx).and_then(|ts| {
                if let topods::TShape::Edge(e) = &**ts {
                    Some(e.clone())
                } else {
                    None
                }
            })
        };
        let Some(ed) = ed else {
            return self.write_edge_curve_by_index_topods(edge_idx);
        };

        let start_pt = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(ed.first.index).and_then(|ts| {
                if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
            }).unwrap_or_default()
        };
        let end_pt = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(ed.last.index).and_then(|ts| {
                if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
            }).unwrap_or_default()
        };

        let seam_is_cylinder = matches!(face_surface.as_ref(), Some(Surface3::Cylinder(_)));
        let seam_is_cone = matches!(face_surface.as_ref(), Some(Surface3::Cone(_)));
        let axis = canonicalize_axis_sign(glam::DVec3::Z);
        let start_proj = start_pt.dot(axis);
        let end_proj = end_pt.dot(axis);
        let canonical_low_to_high = seam_is_cylinder;

        // Get pcurves from edge representations directly
        struct PCurveEntry { face_tshape_idx: usize, curve2d: Curve2d }
        let mut pcurve_entries: Vec<PCurveEntry> = Vec::new();
        for rep in &ed.representations {
            match rep {
                topods::CurveRepresentation::CurveOnSurface { face, pcurve, .. } => {
                    pcurve_entries.push(PCurveEntry { face_tshape_idx: *face, curve2d: pcurve.clone() });
                }
                topods::CurveRepresentation::CurveOnClosedSurface { face, pcurve1, .. } => {
                    pcurve_entries.push(PCurveEntry { face_tshape_idx: *face, curve2d: pcurve1.clone() });
                }
                _ => {}
            }
        }
        let pcurves: Vec<Curve2d> = pcurve_entries.iter().map(|e| e.curve2d.clone()).collect();

        let synthetic_curve2d = if pcurves.is_empty() {
            match face_surface.as_ref() {
                Some(Surface3::Cylinder(cyl)) => {
                    synthesize_cylinder_pcurve_for_edge_topods(self.tbrep(), edge_idx, cyl)
                }
                _ => None,
            }
        } else {
            None
        };

        let basis_curve = match face_surface.clone() {
            Some(Surface3::Sphere(sphere)) => {
                let a = (start_pt - sphere.center).normalize_or_zero();
                let b = (end_pt - sphere.center).normalize_or_zero();
                let mut circle_normal = a.cross(b);
                if circle_normal.length_squared() < 1e-12 {
                    circle_normal = any_perpendicular_dvec3(sphere.axis);
                }
                let circle_normal = circle_normal.normalize_or_zero();
                let placement =
                    self.axis2_from_origin_axis("seam_axis", sphere.center, circle_normal);
                self.circle("seam_circle", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(_cone)) => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Cylinder(_)) => {
                let (origin_pt, delta_vec) = if start_proj <= end_proj {
                    (start_pt, end_pt - start_pt)
                } else {
                    (end_pt, start_pt - end_pt)
                };
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(origin_pt));
                let delta = dvec3_to_array(delta_vec);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Torus(torus)) => {
                let axis = torus.axis.normalize_or_zero();
                if axis.length_squared() < 1e-18 {
                    let origin_id =
                        self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                    let delta = dvec3_to_array(end_pt - start_pt);
                    let magnitude = vector_length(delta).max(1e-9);
                    let dir = self.direction("seam_dir", normalize(delta));
                    let vec = self.vector("seam_vec", dir, magnitude);
                    self.line("seam_line", origin_id, vec)
                } else {
                    let major_seam = pcurves.iter().any(|c2| match c2 {
                        Curve2d::Line(l) => l.direction.x.abs() > l.direction.y.abs(),
                        _ => false,
                    });

                    if major_seam {
                        let radial_vec = start_pt
                            - torus.center
                            - axis * (start_pt - torus.center).dot(axis);
                        let seam_radius = radial_vec.length().max(1e-9);
                        let placement =
                            self.axis2_from_origin_axis("seam_torus_axis", torus.center, axis);
                        self.circle("seam_torus_circle", placement, seam_radius)
                    } else {
                        let mid = (start_pt + end_pt) * 0.5;
                        let radial_raw =
                            mid - torus.center - axis * (mid - torus.center).dot(axis);
                        let radial = if radial_raw.length_squared() < 1e-18 {
                            any_perpendicular_dvec3(axis)
                        } else {
                            radial_raw.normalize_or_zero()
                        };
                        let circle_center = torus.center + radial * torus.major_radius;
                        let circle_normal = axis.cross(radial).normalize_or_zero();
                        if circle_normal.length_squared() < 1e-18 {
                            let origin_id =
                                self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                            let delta = dvec3_to_array(end_pt - start_pt);
                            let magnitude = vector_length(delta).max(1e-9);
                            let dir = self.direction("seam_dir", normalize(delta));
                            let vec = self.vector("seam_vec", dir, magnitude);
                            self.line("seam_line", origin_id, vec)
                        } else {
                            let placement = self.axis2_from_origin_axis(
                                "seam_torus_axis",
                                circle_center,
                                circle_normal,
                            );
                            self.circle(
                                "seam_torus_circle",
                                placement,
                                torus.minor_radius.max(1e-9),
                            )
                        }
                    }
                }
            }
            _ => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
        };

        let final_curve = if !pcurves.is_empty() || synthetic_curve2d.is_some() {
            let mut pcurve_ids = Vec::new();
            let mut periodic_extra_curve2d: Option<Curve2d> = None;
            let first_curve2d = if let Some(pce0) = pcurve_entries.first() {
                Some(pce0.curve2d.clone())
            } else {
                synthetic_curve2d.clone()
            };
            if (seam_is_cylinder || seam_is_cone)
                && (pcurves.len() == 1
                    || (pcurves.is_empty() && synthetic_curve2d.is_some()))
                && let Some(Curve2d::Line(l0)) = first_curve2d
            {
                let eps = 1e-9;
                if l0.direction.x.abs() <= eps && l0.direction.y.abs() > eps {
                    let mut l1 = l0;
                    l1.origin.x = l0.origin.x + 2.0 * std::f64::consts::PI;
                    periodic_extra_curve2d = Some(Curve2d::Line(l1));
                }
            }

            if let Some(curve2d) = synthetic_curve2d
                && let Some(Surface3::Cylinder(cyl)) = face_surface.as_ref()
            {
                let surface_id = self.write_surface(Some(Surface3::Cylinder(*cyl)), None);
                let param_curve_id = self.write_curve2d(Some(curve2d));
                let def_rep = self.definitional_representation(param_curve_id);
                let pcurve_id = self.pcurve(surface_id, def_rep);
                pcurve_ids.push(pcurve_id);
                if let Some(extra_curve2d) = periodic_extra_curve2d.clone() {
                    let extra_param = self.write_curve2d(Some(extra_curve2d));
                    let extra_def = self.definitional_representation(extra_param);
                    let extra_pc = self.pcurve(surface_id, extra_def);
                    pcurve_ids.push(extra_pc);
                }
            }

            for (pc_i, pce) in pcurve_entries.iter().enumerate() {
                let surface_id = self.get_or_write_surface_for_face(pce.face_tshape_idx);
                let mut curve2d = Some(pce.curve2d.clone());
                if seam_is_cylinder
                    && let Some(Curve2d::Line(mut l)) = curve2d
                {
                    let eps = 1e-9;
                    if l.direction.x.abs() <= eps && l.direction.y.abs() > eps {
                        if l.direction.y < 0.0 {
                            l.direction.y = -l.direction.y;
                        }
                        l.origin.y = 0.0;
                        let two_pi = 2.0 * std::f64::consts::PI;
                        let mut u = l.origin.x.rem_euclid(two_pi);
                        if u <= 1e-6 {
                            u = 0.0;
                        }
                        if (two_pi - u).abs() <= 1e-6 {
                            u = two_pi;
                        }
                        l.origin.x = u;
                    }
                    curve2d = Some(Curve2d::Line(l));
                }

                let param_curve_id = self.write_curve2d(curve2d);
                let def_rep = self.definitional_representation(param_curve_id);
                let pcurve_id = self.pcurve(surface_id, def_rep);
                pcurve_ids.push(pcurve_id);

                if pc_i == 0 && let Some(extra_curve2d) = periodic_extra_curve2d.clone() {
                    let extra_param = self.write_curve2d(Some(extra_curve2d));
                    let extra_def = self.definitional_representation(extra_param);
                    let extra_pc = self.pcurve(surface_id, extra_def);
                    pcurve_ids.push(extra_pc);
                }
            }
            if pcurve_ids.len() >= 2 {
                self.seam_curve(basis_curve, &pcurve_ids)
            } else {
                self.surface_curve(basis_curve, &pcurve_ids)
            }
        } else {
            basis_curve
        };

        let v0 = self.vertex_point_by_tshape_idx(ed.first.index);
        let v1 = self.vertex_point_by_tshape_idx(ed.last.index);
        self.edge_curve("seam_edge", v0, v1, final_curve, true)
    }

    fn write_seam_edge_curve(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        face_surface: Option<Surface3>,
    ) -> u64 {
        let Some(edge) = brep.edges.get(edge_idx) else {
            return self.write_edge_curve_by_index_topods(edge_idx);
        };
        let start_pt = brep
            .vertices
            .get(edge.start)
            .map(|v| v.point)
            .unwrap_or(glam::DVec3::ZERO);
        let end_pt = brep
            .vertices
            .get(edge.end)
            .map(|v| v.point)
            .unwrap_or(glam::DVec3::ZERO);

        let seam_is_cylinder = matches!(face_surface.as_ref(), Some(Surface3::Cylinder(_)));
        let seam_is_cone = matches!(face_surface.as_ref(), Some(Surface3::Cone(_)));
        let axis = canonicalize_axis_sign(glam::DVec3::Z);
        let start_proj = start_pt.dot(axis);
        let end_proj = end_pt.dot(axis);
        let canonical_low_to_high = seam_is_cylinder;

        let (edge_start_idx, edge_end_idx) = if canonical_low_to_high && start_proj > end_proj {
            (edge.end, edge.start)
        } else {
            (edge.start, edge.end)
        };
        let v0 = self.vertex_point_by_index(brep, edge_start_idx);
        let v1 = self.vertex_point_by_index(brep, edge_end_idx);

        let pcurves = brep
            .geom
            .edge_pcurves
            .get(edge_idx)
            .cloned()
            .unwrap_or_default();

        let synthetic_curve2d = if pcurves.is_empty() {
            match face_surface.as_ref() {
                Some(Surface3::Cylinder(cyl)) => synthesize_cylinder_pcurve_for_edge(brep, edge_idx, cyl),
                _ => None,
            }
        } else {
            None
        };

        let basis_curve = match face_surface.clone() {
            Some(Surface3::Sphere(sphere)) => {
                let a = (start_pt - sphere.center).normalize_or_zero();
                let b = (end_pt - sphere.center).normalize_or_zero();
                let mut circle_normal = a.cross(b);
                if circle_normal.length_squared() < 1e-12 {
                    circle_normal = any_perpendicular_dvec3(sphere.axis);
                }
                let circle_normal = circle_normal.normalize_or_zero();
                let placement =
                    self.axis2_from_origin_axis("seam_axis", sphere.center, circle_normal);
                self.circle("seam_circle", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(_cone)) => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Cylinder(_)) => {
                let (origin_pt, delta_vec) = if start_proj <= end_proj {
                    (start_pt, end_pt - start_pt)
                } else {
                    (end_pt, start_pt - end_pt)
                };
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(origin_pt));
                let delta = dvec3_to_array(delta_vec);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Torus(torus)) => {
                let axis = torus.axis.normalize_or_zero();
                if axis.length_squared() < 1e-18 {
                    let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                    let delta = dvec3_to_array(end_pt - start_pt);
                    let magnitude = vector_length(delta).max(1e-9);
                    let dir = self.direction("seam_dir", normalize(delta));
                    let vec = self.vector("seam_vec", dir, magnitude);
                    self.line("seam_line", origin_id, vec)
                } else {
                    let major_seam = pcurves.iter().any(|pc| {
                        brep
                            .geom
                            .curve2ds
                            .get(pc.curve2d_idx)
                            .cloned()
                            .and_then(|c2| match c2 {
                                Curve2d::Line(l) => {
                                    Some(l.direction.x.abs() > l.direction.y.abs())
                                }
                                _ => None,
                            })
                            .unwrap_or(false)
                    });

                    if major_seam {
                        let radial_vec =
                            start_pt - torus.center - axis * (start_pt - torus.center).dot(axis);
                        let seam_radius = radial_vec.length().max(1e-9);
                        let placement =
                            self.axis2_from_origin_axis("seam_torus_axis", torus.center, axis);
                        self.circle("seam_torus_circle", placement, seam_radius)
                    } else {
                        let mid = (start_pt + end_pt) * 0.5;
                        let radial_raw =
                            mid - torus.center - axis * (mid - torus.center).dot(axis);
                        let radial = if radial_raw.length_squared() < 1e-18 {
                            any_perpendicular_dvec3(axis)
                        } else {
                            radial_raw.normalize_or_zero()
                        };
                        let circle_center = torus.center + radial * torus.major_radius;
                        let circle_normal = axis.cross(radial).normalize_or_zero();
                        if circle_normal.length_squared() < 1e-18 {
                            let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                            let delta = dvec3_to_array(end_pt - start_pt);
                            let magnitude = vector_length(delta).max(1e-9);
                            let dir = self.direction("seam_dir", normalize(delta));
                            let vec = self.vector("seam_vec", dir, magnitude);
                            self.line("seam_line", origin_id, vec)
                        } else {
                            let placement = self.axis2_from_origin_axis(
                                "seam_torus_axis",
                                circle_center,
                                circle_normal,
                            );
                            self.circle("seam_torus_circle", placement, torus.minor_radius.max(1e-9))
                        }
                    }
                }
            }
            _ => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
        };

        let final_curve = if !pcurves.is_empty() || synthetic_curve2d.is_some() {
            let mut pcurve_ids = Vec::new();
            let mut periodic_extra_curve2d: Option<Curve2d> = None;
            let first_curve2d = if let Some(pc0) = pcurves.first() {
                brep.geom.curve2ds.get(pc0.curve2d_idx).cloned()
            } else {
                synthetic_curve2d.clone()
            };
            if (seam_is_cylinder || seam_is_cone)
                && (pcurves.len() == 1 || (pcurves.is_empty() && synthetic_curve2d.is_some()))
                && let Some(Curve2d::Line(l0)) = first_curve2d
                {
                    let eps = 1e-9;
                    if l0.direction.x.abs() <= eps && l0.direction.y.abs() > eps {
                        let mut l1 = l0;
                        l1.origin.x = l0.origin.x + 2.0 * std::f64::consts::PI;
                        periodic_extra_curve2d = Some(Curve2d::Line(l1));
                    }
                }

            if let Some(curve2d) = synthetic_curve2d
                && let Some(Surface3::Cylinder(cyl)) = face_surface.as_ref()
            {
                let surface_id = self.write_surface(Some(Surface3::Cylinder(*cyl)), None);
                let param_curve_id = self.write_curve2d(Some(curve2d));
                let def_rep = self.definitional_representation(param_curve_id);
                let pcurve_id = self.pcurve(surface_id, def_rep);
                pcurve_ids.push(pcurve_id);
                if let Some(extra_curve2d) = periodic_extra_curve2d.clone() {
                    let extra_param = self.write_curve2d(Some(extra_curve2d));
                    let extra_def = self.definitional_representation(extra_param);
                    let extra_pc = self.pcurve(surface_id, extra_def);
                    pcurve_ids.push(extra_pc);
                }
            }

            for (pc_i, pc) in pcurves.iter().enumerate() {
                let surface_id = self.get_or_write_surface_id(brep, pc.surface_idx);
                let mut curve2d = brep.geom.curve2ds.get(pc.curve2d_idx).cloned();
                if seam_is_cylinder
                    && let Some(Curve2d::Line(mut l)) = curve2d
                {
                    let eps = 1e-9;
                    if l.direction.x.abs() <= eps && l.direction.y.abs() > eps {
                        if l.direction.y < 0.0 {
                            l.direction.y = -l.direction.y;
                        }
                        l.origin.y = 0.0;
                        let two_pi = 2.0 * std::f64::consts::PI;
                        let mut u = l.origin.x.rem_euclid(two_pi);
                        if u <= 1e-6 {
                            u = 0.0;
                        }
                        if (two_pi - u).abs() <= 1e-6 {
                            u = two_pi;
                        }
                        l.origin.x = u;
                    }
                    curve2d = Some(Curve2d::Line(l));
                }

                let param_curve_id = self.write_curve2d(curve2d);
                let def_rep = self.definitional_representation(param_curve_id);
                let pcurve_id = self.pcurve(surface_id, def_rep);
                pcurve_ids.push(pcurve_id);

                if pc_i == 0 && let Some(extra_curve2d) = periodic_extra_curve2d.clone() {
                    let extra_param = self.write_curve2d(Some(extra_curve2d));
                    let extra_def = self.definitional_representation(extra_param);
                    let extra_pc = self.pcurve(surface_id, extra_def);
                    pcurve_ids.push(extra_pc);
                }
            }
            if pcurve_ids.len() >= 2 {
                self.seam_curve(basis_curve, &pcurve_ids)
            } else {
                self.surface_curve(basis_curve, &pcurve_ids)
            }
        } else {
            basis_curve
        };

        self.edge_curve("seam_edge", v0, v1, final_curve, true)
    }

    fn write_surface(
        &mut self,
        face_surface: Option<Surface3>,
        fallback_placement: Option<u64>,
    ) -> u64 {
        match face_surface {
            Some(Surface3::Plane(plane)) => {
                let plane_normal = canonicalize_axis_sign(plane.normal);
                let placement =
                    self.axis2_from_origin_axis("plane_axis", plane.origin, plane_normal);
                self.plane("face_plane", placement)
            }
            Some(Surface3::Cylinder(cyl)) => {
                let placement = self.axis2_from_origin_axis("cyl_axis", cyl.origin, cyl.axis);
                self.cylindrical_surface("face_cylinder", placement, cyl.radius.max(1e-9))
            }
            Some(Surface3::Sphere(sphere)) => {
                let placement =
                    self.axis2_from_origin_axis("sphere_axis", sphere.center, sphere.axis);
                self.spherical_surface("face_sphere", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(cone)) => {
                let placement = self.axis2_from_origin_axis("cone_axis", cone.apex, cone.axis);
                // STEP semi-angle is in radians per AP214/AP242 standard.
                self.conical_surface(
                    "face_cone",
                    placement,
                    cone.radius,
                    cone.half_angle_rad,
                )
            }
            Some(Surface3::Torus(torus)) => {
                let placement = self.axis2_from_origin_axis("torus_axis", torus.center, torus.axis);
                self.toroidal_surface(
                    "face_torus",
                    placement,
                    torus.major_radius.max(1e-9),
                    torus.minor_radius.max(1e-9),
                )
            }
            None => {
                let placement = if let Some(id) = fallback_placement {
                    id
                } else {
                    let origin = self.cartesian_point("face_origin", [0.0, 0.0, 0.0]);
                    let axis = self.direction("face_normal", [0.0, 0.0, 1.0]);
                    let ref_dir = self.direction("face_ref", [1.0, 0.0, 0.0]);
                    self.axis2_placement_3d("face_axis", origin, axis, ref_dir)
                };
                self.plane("face_plane", placement)
            }
            Some(Surface3::BSpline(bs)) => {
                // 鈿狅笍 淇濇寔 BSpline 琛ㄩ潰绫诲瀷涓嶅彉,涓嶆彁鍗囦负 PLANE锟?
                //    OCCT 鍙傦拷?STEP 鏂囦欢淇濈暀鍘熷琛ㄩ潰绫诲瀷(锟?NURBS 杞崲鍚庣殑 BSpline 锟?锟?
                self.write_bspline_surface(&bs.clone())
            }
            Some(Surface3::Ellipsoid(ellipsoid)) => {
                let bs = surface_to_bspline(&Surface3::Ellipsoid(ellipsoid), 9, 9);
                let name = format!(
                    "RCAD_ELLIPSOID;c={},{},{};a={},{},{};d={},{},{};r={},{},{}",
                    ellipsoid.center.x,
                    ellipsoid.center.y,
                    ellipsoid.center.z,
                    ellipsoid.axis.x,
                    ellipsoid.axis.y,
                    ellipsoid.axis.z,
                    ellipsoid.ref_dir.x,
                    ellipsoid.ref_dir.y,
                    ellipsoid.ref_dir.z,
                    ellipsoid.radius_x,
                    ellipsoid.radius_y,
                    ellipsoid.radius_z,
                );
                self.write_bspline_surface_named(&name, &bs)
            }
            Some(Surface3::Helicoid(helicoid)) => {
                let bs = surface_to_bspline(&Surface3::Helicoid(helicoid), 9, 9);
                let name = format!(
                    "RCAD_HELICOID;o={},{},{};a={},{},{};d={},{},{};p={}",
                    helicoid.origin.x,
                    helicoid.origin.y,
                    helicoid.origin.z,
                    helicoid.axis.x,
                    helicoid.axis.y,
                    helicoid.axis.z,
                    helicoid.ref_dir.x,
                    helicoid.ref_dir.y,
                    helicoid.ref_dir.z,
                    helicoid.pitch,
                );
                self.write_bspline_surface_named(&name, &bs)
            }
            Some(Surface3::Revolution(rev)) => {
                // Write as SURFACE_OF_REVOLUTION referencing the profile curve
                // and axis placement, matching OCCT STEP output.
                let domain = rev.profile.default_domain();
                let p_start = if domain[0].is_finite() {
                    rev.profile.point_at(domain[0])
                } else {
                    glam::DVec3::ZERO
                };
                let p_end = if domain[1].is_finite() {
                    rev.profile.point_at(domain[1])
                } else {
                    glam::DVec3::new(0.0, 0.0, 1.0)
                };
                let start_pt = dvec3_to_array(p_start);
                let end_pt = dvec3_to_array(p_end);
                let curve_id = self.write_curve3_entity(&rev.profile, start_pt, end_pt);
                let origin = self.cartesian_point("rev_origin", dvec3_to_array(rev.axis_origin));
                let axis = self.direction("rev_axis", normalize(dvec3_to_array(rev.axis_dir)));
                let axis_placement = self.axis1_placement("rev_axis_placement", origin, axis);
                self.surface_of_revolution("face_revolution", curve_id, axis_placement)
            }
            Some(Surface3::LinearExtrusion(ext)) => {
                // Write as SURFACE_OF_LINEAR_EXTRUSION referencing the profile
                // curve and extrusion vector, matching OCCT STEP output.
                let domain = ext.profile.default_domain();
                let p_start = if domain[0].is_finite() {
                    ext.profile.point_at(domain[0])
                } else {
                    glam::DVec3::ZERO
                };
                let p_end = if domain[1].is_finite() {
                    ext.profile.point_at(domain[1])
                } else {
                    glam::DVec3::new(0.0, 0.0, 1.0)
                };
                let start_pt = dvec3_to_array(p_start);
                let end_pt = dvec3_to_array(p_end);
                let curve_id = self.write_curve3_entity(&ext.profile, start_pt, end_pt);
                let dir = normalize(dvec3_to_array(ext.direction));
                let dir_id = self.direction("ext_dir", dir);
                let vec_id = self.vector("ext_vec", dir_id, 1.0);
                self.surface_of_linear_extrusion("face_extrusion", curve_id, vec_id)
            }
            Some(Surface3::Offset(offset)) => {
                // Write as OFFSET_SURFACE referencing the basis surface and
                // offset distance, matching OCCT STEP output.
                let basis_id = self.write_surface(Some(*offset.basis), None);
                self.offset_surface("face_offset", basis_id, offset.offset_distance)
            }
            Some(surface @ Surface3::Pipe(_))
            | Some(surface @ Surface3::Ruled(_))
            | Some(surface @ Surface3::Coons(_))
            | Some(surface @ Surface3::TriBezier(_))
            | Some(surface @ Surface3::Bezier(_)) => {
                // Export higher-level surfaces through a sampled NURBS fallback
                // instead of collapsing them to a plane.
                let bs = surface_to_bspline(&surface, 9, 9);
                self.write_bspline_surface(&bs)
            }
            Some(Surface3::Trimmed(ts)) => {
                // Write the underlying basis surface; trim bounds are implied by the
                // face topology, so we don't need a separate RECTANGULAR_TRIMMED_SURFACE entity.
                self.write_surface(Some(*ts.basis), fallback_placement)
            }
        }
    }

    fn write_bspline_surface(&mut self, bs: &rcad_kernel::geom::BSplineSurface) -> u64 {
        self.write_bspline_surface_named("bspline_surf", bs)
    }

    fn write_bspline_surface_named(
        &mut self,
        name: &str,
        bs: &rcad_kernel::geom::BSplineSurface,
    ) -> u64 {
        let n_u = bs.control_points.len();
        let n_v = bs.control_points.first().map(|r| r.len()).unwrap_or(0);
        if n_u == 0 || n_v == 0 {
            let origin = self.cartesian_point("bs_origin", [0.0, 0.0, 0.0]);
            let ax = self.direction("bs_norm", [0.0, 0.0, 1.0]);
            let rd = self.direction("bs_ref", [1.0, 0.0, 0.0]);
            let pl = self.axis2_placement_3d("bs_pl", origin, ax, rd);
            return self.plane("bs_fallback", pl);
        }
        // Build cp_grid[v][u] for STEP (STEP stores [v][u], we store [u][v])
        let cp_grid: Vec<Vec<u64>> = (0..n_v)
            .map(|vi| {
                (0..n_u)
                    .map(|ui| {
                        self.cartesian_point("bs_cp", dvec3_to_array(bs.control_points[ui][vi]))
                    })
                    .collect()
            })
            .collect();
        let (mults_u, knots_u) = compress_knot_vector(&bs.knots_u);
        let (mults_v, knots_v) = compress_knot_vector(&bs.knots_v);
        self.b_spline_surface_with_knots(
            name,
            bs.degree_u,
            bs.degree_v,
            &cp_grid,
            &mults_u,
            &knots_u,
            &mults_v,
            &knots_v,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn b_spline_surface_with_knots(
        &mut self,
        name: &str,
        degree_u: usize,
        degree_v: usize,
        cp_grid: &[Vec<u64>],
        mults_u: &[usize],
        knots_u: &[f64],
        mults_v: &[usize],
        knots_v: &[f64],
    ) -> u64 {
        let rows: Vec<String> = cp_grid
            .iter()
            .map(|row| {
                let refs = row
                    .iter()
                    .map(|&r| format!("#{}", r))
                    .collect::<Vec<_>>()
                    .join(",");
                format!("({})", refs)
            })
            .collect();
        let cp_str = format!("({})", rows.join(","));
        let mu_str: String = format!(
            "({})",
            mults_u
                .iter()
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        let mv_str: String = format!(
            "({})",
            mults_v
                .iter()
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        let ku_str: String = format!(
            "({})",
            knots_u
                .iter()
                .map(|k| format!("{:.9}", k))
                .collect::<Vec<_>>()
                .join(",")
        );
        let kv_str: String = format!(
            "({})",
            knots_v
                .iter()
                .map(|k| format!("{:.9}", k))
                .collect::<Vec<_>>()
                .join(",")
        );
        self.push(format!(
            "B_SPLINE_SURFACE_WITH_KNOTS('{}',{},{},{},.UNSPECIFIED.,.F.,.F.,.F.,{},{},{},{},.UNSPECIFIED.)",
            name, degree_u, degree_v, cp_str, mu_str, mv_str, ku_str, kv_str
        ))
    }

    fn axis2_from_origin_axis(
        &mut self,
        name: &str,
        origin: glam::DVec3,
        axis: glam::DVec3,
    ) -> u64 {
        let origin_id = self.cartesian_point("surface_origin", dvec3_to_array(origin));
        let axis_arr = normalize(dvec3_to_array(axis));
        let axis_id = self.direction("surface_axis", axis_arr);
        let ref_dir = if axis_arr[2] > 0.999999 {
            [1.0, 0.0, 0.0]
        } else {
            orthogonal_dir(axis_arr)
        };
        let ref_id = self.direction("surface_ref", ref_dir);
        self.axis2_placement_3d(name, origin_id, axis_id, ref_id)
    }

    fn axis2_from_origin_axis_ref(
        &mut self,
        name: &str,
        origin: glam::DVec3,
        axis: glam::DVec3,
        ref_dir: glam::DVec3,
    ) -> u64 {
        let origin_id = self.cartesian_point("curve_origin", dvec3_to_array(origin));
        let axis_arr = normalize(dvec3_to_array(axis));
        let axis_id = self.direction("curve_axis", axis_arr);
        let ref_arr = normalize(project_to_plane(dvec3_to_array(ref_dir), axis_arr));
        let ref_id = self.direction("curve_ref", ref_arr);
        self.axis2_placement_3d(name, origin_id, axis_id, ref_id)
    }

    /// Topods-native version of write_edge_curve_by_index (WIP).
    /// Not yet wired into the call chain — all callers still use the old flat version.
    fn write_edge_curve_by_index_topods(&mut self, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_curve_ids.get(&edge_idx) {
            return *existing;
        }
        let ed = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(edge_idx).and_then(|ts| {
                if let topods::TShape::Edge(e) = &**ts { Some(e.clone()) } else { None }
            })
        };
        let Some(ed) = ed else {
            let p0 = self.cartesian_point("edge_p0", [0.0, 0.0, 0.0]);
            let p1 = self.cartesian_point("edge_p1", [0.0, 0.0, 0.0]);
            let ve0 = self.vertex_point("v0", p0);
            let ve1 = self.vertex_point("v1", p1);
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            let basis = self.line("edge_line", origin, vec);
            let ec = self.edge_curve("edge", ve0, ve1, basis, true);
            self.edge_curve_ids.insert(edge_idx, ec);
            return ec;
        };
        let v0 = self.vertex_point_by_tshape_idx(ed.first.index);
        let v1 = self.vertex_point_by_tshape_idx(ed.last.index);
        let curve = ed.curve.clone();
        let start_pt = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(ed.first.index).and_then(|ts| {
                if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
            }).unwrap_or_default()
        };
        let end_pt = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(ed.last.index).and_then(|ts| {
                if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
            }).unwrap_or_default()
        };
        let basis_curve = self.write_edge_curve_geometry_by_index_topods(edge_idx);
        let mut same_sense = true;
        if let Some((center, axis, major_dir)) = match curve {
            Some(Curve3::Circle(c)) => Some((c.center, c.normal, glam::DVec3::ZERO)),
            Some(Curve3::Ellipse(e)) => Some((e.center, e.normal, e.major_dir)),
            _ => None,
        } {
            let canon_axis = canonicalize_axis_sign(axis);
            let ref_dir = if major_dir.length_squared() > 1e-30 {
                let proj = major_dir - major_dir.dot(canon_axis) * canon_axis;
                let len = proj.length();
                if len > 1e-15 { proj / len } else { canon_axis.cross(glam::DVec3::Y).normalize_or_zero() }
            } else if canon_axis.z.abs() > 0.999999 { glam::DVec3::X }
            else {
                let helper = if canon_axis.y.abs() < 0.9 { glam::DVec3::Y } else { glam::DVec3::X };
                canon_axis.cross(helper).normalize()
            };
            let perp = canon_axis.cross(ref_dir);
            let d_start = start_pt - center;
            let d_end = end_pt - center;
            let theta_start = f64::atan2(d_start.dot(perp), d_start.dot(ref_dir));
            let theta_end = f64::atan2(d_end.dot(perp), d_end.dot(ref_dir));
            let theta_start = theta_start.rem_euclid(std::f64::consts::TAU);
            let theta_end = theta_end.rem_euclid(std::f64::consts::TAU);
            let forward = if theta_end >= theta_start {
                theta_end - theta_start
            } else {
                theta_end + std::f64::consts::TAU - theta_start
            };
            if forward > std::f64::consts::PI { same_sense = false; }
        }
        let edge_curve = self.edge_curve("edge", v0, v1, basis_curve, same_sense);
        self.edge_curve_ids.insert(edge_idx, edge_curve);
        edge_curve
    }

    fn write_edge_curve_by_index(&mut self, brep: &BRep, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_curve_ids.get(&edge_idx) {
            return *existing;
        }

        let Some(edge) = brep.edges.get(edge_idx) else {
            let p0 = self.cartesian_point("edge_p0", [0.0, 0.0, 0.0]);
            let p1 = self.cartesian_point("edge_p1", [0.0, 0.0, 0.0]);
            let v0 = self.vertex_point("v0", p0);
            let v1 = self.vertex_point("v1", p1);
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            let basis = self.line("edge_line", origin, vec);
            return self.edge_curve("edge", v0, v1, basis, true);
        };

        let _start_point = brep
            .vertices
            .get(edge.start)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let _end_point = brep
            .vertices
            .get(edge.end)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let v0 = self.vertex_point_by_index(brep, edge.start);
        let v1 = self.vertex_point_by_index(brep, edge.end);
        let basis_curve = self.write_edge_curve_geometry_by_index(brep, edge_idx);

        // Determine `same_sense` for periodic edge curves (CIRCLE, ELLIPSE).
        // The STEP CIRCLE / ELLIPSE entity uses `axis2_from_origin_axis` which
        // picks a specific `ref_dir` that may differ from the BRep curve's
        // internal parameterization.  We compute the counterclockwise angle
        // from start鈫抏nd on the curve's frame; if that angle exceeds 蟺 we
        // set `same_sense = false` so the STEP reader takes the complement
        // (shorter) arc, avoiding ambiguous parametrisation.
        let mut same_sense = true;
        let curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .copied()
            .flatten()
            .and_then(|ci| brep.geom.curves.get(ci));
        if let Some((center, axis, major_dir)) = match curve {
            Some(Curve3::Circle(c)) => Some((c.center, c.normal, glam::DVec3::ZERO)),
            Some(Curve3::Ellipse(e)) => Some((e.center, e.normal, e.major_dir)),
            _ => None,
        } {
            let canon_axis = canonicalize_axis_sign(axis);
            // ref_dir matches axis2_from_origin_axis (CIRCLE) /
            // axis2_from_origin_axis_ref with major_dir (ELLIPSE).
            let ref_dir = if major_dir.length_squared() > 1e-30 {
                // ELLIPSE: use the stored major_dir projected onto the plane.
                let proj = major_dir - major_dir.dot(canon_axis) * canon_axis;
                let len = proj.length();
                if len > 1e-15 { proj / len } else { canon_axis.cross(glam::DVec3::Y).normalize_or_zero() }
            } else if canon_axis.z.abs() > 0.999999 {
                glam::DVec3::X
            } else {
                let helper = if canon_axis.y.abs() < 0.9 {
                    glam::DVec3::Y
                } else {
                    glam::DVec3::X
                };
                canon_axis.cross(helper).normalize()
            };
            let perp = canon_axis.cross(ref_dir);
            let start_pt = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or_default();
            let end_pt = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or_default();
            let d_start = start_pt - center;
            let d_end = end_pt - center;
            let theta_start = f64::atan2(d_start.dot(perp), d_start.dot(ref_dir));
            let theta_end = f64::atan2(d_end.dot(perp), d_end.dot(ref_dir));
            let theta_start = theta_start.rem_euclid(std::f64::consts::TAU);
            let theta_end = theta_end.rem_euclid(std::f64::consts::TAU);
            let forward = if theta_end >= theta_start {
                theta_end - theta_start
            } else {
                theta_end + std::f64::consts::TAU - theta_start
            };
            if forward > std::f64::consts::PI {
                same_sense = false;
            }
        }
        let edge_curve = self.edge_curve("edge", v0, v1, basis_curve, same_sense);
        self.edge_curve_ids.insert(edge_idx, edge_curve);
        edge_curve
    }

    /// Write surface entity for a face tshape index (topods version of get_or_write_surface_id).
    fn get_or_write_surface_for_face(&mut self, face_idx: usize) -> u64 {
        if let Some(existing) = self.surface_ids.get(&face_idx) {
            return *existing;
        }
        let surface = self.tbrep().tshapes.get(face_idx).and_then(|ts| {
            if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
        });
        let sid = self.write_surface(surface, None);
        self.surface_ids.insert(face_idx, sid);
        sid
    }

    /// Topods variant: same logic as write_edge_curve_geometry_by_index but reads
    /// geometry data from self.tbrep() (topods::BRep) instead of flat BRep.
    fn write_edge_curve_geometry_by_index_topods(&mut self, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_geometry_ids.get(&edge_idx) {
            return *existing;
        }

        let ed_opt = {
            let tbrep = self.tbrep();
            tbrep.tshapes.get(edge_idx).and_then(|ts| {
                if let topods::TShape::Edge(e) = &**ts { Some(e.clone()) } else { None }
            })
        };

        let curve_id = if ed_opt.is_none() {
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            self.line("edge_line", origin, vec)
        } else {
            let ed = ed_opt.as_ref().unwrap();

            let start_point = {
                let tbrep = self.tbrep();
                tbrep.tshapes.get(ed.first.index).and_then(|ts| {
                    if let topods::TShape::Vertex(v) = &**ts { Some(dvec3_to_array(v.point)) } else { None }
                }).unwrap_or([0.0, 0.0, 0.0])
            };
            let end_point = {
                let tbrep = self.tbrep();
                tbrep.tshapes.get(ed.last.index).and_then(|ts| {
                    if let topods::TShape::Vertex(v) = &**ts { Some(dvec3_to_array(v.point)) } else { None }
                }).unwrap_or([0.0, 0.0, 0.0])
            };

            let source_curve = ed.curve.clone();

            // Build pcurves from the edge CurveRepresentations.
            // surface_idx in each PCurve is the face tshape index.
            // curve2d_idx indexes into rep_curve2ds.
            let mut pcurves: Vec<flat::PCurve> = Vec::new();
            let mut rep_curve2ds: Vec<Curve2d> = Vec::new();
            for rep in &ed.representations {
                match rep {
                    topods::CurveRepresentation::CurveOnSurface { face, pcurve, .. } => {
                        let ci = rep_curve2ds.len();
                        rep_curve2ds.push(pcurve.clone());
                        pcurves.push(flat::PCurve { surface_idx: *face, curve2d_idx: ci });
                    }
                    topods::CurveRepresentation::CurveOnClosedSurface { face, pcurve1, pcurve2, .. } => {
                        let ci1 = rep_curve2ds.len();
                        rep_curve2ds.push(pcurve1.clone());
                        pcurves.push(flat::PCurve { surface_idx: *face, curve2d_idx: ci1 });
                        let ci2 = rep_curve2ds.len();
                        rep_curve2ds.push(pcurve2.clone());
                        pcurves.push(flat::PCurve { surface_idx: *face, curve2d_idx: ci2 });
                    }
                    _ => {}
                }
            }

            // Inline write_basis_curve_for_edge: write 3d curve from source_curve.
            let mut basis_curve_3d = match source_curve.clone() {
                Some(c) => self.write_curve3_entity(&c, start_point, end_point),
                None => {
                    let p0 = self.cartesian_point("edge_origin", start_point);
                    let delta = [
                        end_point[0] - start_point[0],
                        end_point[1] - start_point[1],
                        end_point[2] - start_point[2],
                    ];
                    let magnitude = vector_length(delta).max(1e-9);
                    let direction = self.direction("edge_dir", normalize(delta));
                    let vector = self.vector("edge_vec", direction, magnitude);
                    self.line("edge_line", p0, vector)
                }
            };

            // Helper: look up surface from a face tshape index (returns owned to avoid lifetime issues).
            let get_surface = |tbrep: &topods::BRep, face_idx: usize| -> Option<Surface3> {
                tbrep.tshapes.get(face_idx).and_then(|ts| {
                    if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
                })
            };

            let is_single_plane_closed_circle = !self.strict_plane_closed_ellipse_done
                && ed.first.index == ed.last.index
                && pcurves.len() == 1
                && matches!(
                    get_surface(self.tbrep(), pcurves[0].surface_idx),
                    Some(Surface3::Plane(_))
                )
                && matches!(source_curve.as_ref(), Some(Curve3::Circle(_)));

            if is_single_plane_closed_circle
                && let Some(Curve3::Circle(circle)) = source_curve.as_ref()
            {
                let placement =
                    self.axis2_from_origin_axis("edge_plane_closed_axis", circle.center, circle.normal);
                let major = circle.radius.max(1e-9);
                let minor = (circle.radius * (1.0 - 1e-9)).max(1e-9);
                basis_curve_3d = self.ellipse("edge_plane_closed_ellipse", placement, major, minor);
                self.strict_plane_closed_ellipse_done = true;
            }

            let torus_surface = if !pcurves.is_empty() {
                let mut torus = None;
                let mut all_torus = true;
                for pc in &pcurves {
                    match get_surface(self.tbrep(), pc.surface_idx) {
                        Some(Surface3::Torus(t)) => {
                            if torus.is_none() {
                                torus = Some(t);
                            }
                        }
                        _ => {
                            all_torus = false;
                            break;
                        }
                    }
                }
                if all_torus { torus } else { None }
            } else {
                None
            };

            if let Some(torus) = torus_surface {
                let mut major_seam = false;
                if let Some(pc0) = pcurves.first()
                    && let Some(Curve2d::Line(l)) = rep_curve2ds.get(pc0.curve2d_idx).cloned()
                {
                    major_seam = l.direction.x.abs() > l.direction.y.abs();
                }

                let start_pt = {
                    let tbrep = self.tbrep();
                    tbrep.tshapes.get(ed.first.index).and_then(|ts| {
                        if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
                    }).unwrap_or(glam::DVec3::ZERO)
                };
                let end_pt = {
                    let tbrep = self.tbrep();
                    tbrep.tshapes.get(ed.last.index).and_then(|ts| {
                        if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
                    }).unwrap_or(glam::DVec3::ZERO)
                };
                let axis = torus.axis.normalize_or_zero();
                if axis.length_squared() >= 1e-18 {
                    if major_seam {
                        let radial_vec = start_pt - torus.center - axis * (start_pt - torus.center).dot(axis);
                        let seam_radius = radial_vec.length().max(1e-9);
                        let placement =
                            self.axis2_from_origin_axis("edge_torus_major_axis", torus.center, axis);
                        basis_curve_3d =
                            self.circle("edge_torus_major_circle", placement, seam_radius);
                    } else {
                        let mid = (start_pt + end_pt) * 0.5;
                        let radial_raw = mid - torus.center - axis * (mid - torus.center).dot(axis);
                        let radial = if radial_raw.length_squared() < 1e-18 {
                            any_perpendicular_dvec3(axis)
                        } else {
                            radial_raw.normalize_or_zero()
                        };
                        let circle_center = torus.center + radial * torus.major_radius;
                        let circle_normal = axis.cross(radial).normalize_or_zero();
                        if circle_normal.length_squared() >= 1e-18 {
                            let placement = self.axis2_from_origin_axis(
                                "edge_torus_seam_axis",
                                circle_center,
                                circle_normal,
                            );
                            basis_curve_3d = self.circle(
                                "edge_torus_seam_circle",
                                placement,
                                torus.minor_radius.max(1e-9),
                            );
                        }
                    }
                }
            }

            if !pcurves.is_empty() {
                let mut pcurve_ids = Vec::new();
                let mut periodic_extra_curve2d: Option<Curve2d> = None;
                let mut periodic_line_dup_for_seam = false;
                if pcurves.len() == 1
                    && let Some(pc0) = pcurves.first() {
                        let is_periodic_u_surface = matches!(
                            get_surface(self.tbrep(), pc0.surface_idx),
                            Some(Surface3::Torus(_))
                                | Some(Surface3::Cylinder(_))
                                | Some(Surface3::Cone(_))
                        );
                        if is_periodic_u_surface
                            && let Some(Curve2d::Line(l0)) = rep_curve2ds.get(pc0.curve2d_idx).cloned() {
                                let eps = 1e-9;
                                let mut l1 = l0;
                                let mut duplicated = false;
                                if l0.direction.x.abs() <= eps && l0.direction.y.abs() > eps {
                                    l1.origin.x = l0.origin.x + 2.0 * std::f64::consts::PI;
                                    duplicated = true;
                                }
                                if duplicated {
                                    periodic_extra_curve2d = Some(Curve2d::Line(l1));
                                    periodic_line_dup_for_seam = true;
                                }
                            }
                    }
                let mut first_plane: Option<rcad_kernel::geom::Plane> = None;
                let mut first_is_line = false;
                for (pc_i, pc) in pcurves.iter().enumerate() {
                    let surface_id = self.get_or_write_surface_for_face(pc.surface_idx);
                    // Get curve2d from rep_curve2ds using pc.curve2d_idx (not from flat GeomStore)
                    let curve2d = rep_curve2ds.get(pc.curve2d_idx).cloned();

                    // curve2d_to_match: use for matching patterns that check the type
                    let mut curve2d = curve2d;
                    let is_cylinder = matches!(
                        get_surface(self.tbrep(), pc.surface_idx),
                        Some(Surface3::Cylinder(_))
                    );

                    if is_cylinder
                        && let Some(Curve2d::Line(mut l)) = curve2d
                    {
                        let eps = 1e-9;
                        if l.direction.y.abs() <= eps && l.direction.x.abs() > eps {
                            l.direction.x = l.direction.x.abs();
                            l.direction.y = 0.0;
                            let two_pi = 2.0 * std::f64::consts::PI;
                            let u = l.origin.x.rem_euclid(two_pi);
                            if u <= 1e-6 || (two_pi - u) <= 1e-6 {
                                l.origin.x = 0.0;
                            } else {
                                l.origin.x = u;
                            }
                        }
                        curve2d = Some(Curve2d::Line(l));
                    }

                    if pcurves.len() == 1 {
                        first_plane = get_surface(self.tbrep(), pc.surface_idx)
                            .and_then(|s| match s {
                                Surface3::Plane(p) => Some(p),
                                _ => None,
                            });
                        first_is_line = matches!(curve2d, Some(Curve2d::Line(_)) | None);
                    }

                    let param_curve_id = self.write_curve2d(curve2d);
                    let def_rep = self.definitional_representation(param_curve_id);
                    let pcurve_id = self.pcurve(surface_id, def_rep);
                    pcurve_ids.push(pcurve_id);

                    if pc_i == 0
                        && let Some(extra_curve2d) = periodic_extra_curve2d.clone()
                    {
                        let extra_param = self.write_curve2d(Some(extra_curve2d));
                        let extra_def = self.definitional_representation(extra_param);
                        let extra_pc = self.pcurve(surface_id, extra_def);
                        pcurve_ids.push(extra_pc);
                    }
                }

                if pcurve_ids.len() == 1
                    && first_is_line
                {
                    let mut extra_surface: Option<(Option<usize>, rcad_kernel::geom::Plane)> = None;
                    if let Some(base_plane) = first_plane
                        && let Some((peer_surface_idx, extra_plane)) =
                            find_peer_plane_surface_for_line_edge_topods(self.tbrep(), edge_idx, &base_plane)
                    {
                        extra_surface = Some((peer_surface_idx, extra_plane));
                    }

                    if extra_surface.is_none()
                        && let Some((surface_idx, plane)) = find_plane_surface_for_edge_topods(self.tbrep(), edge_idx)
                    {
                        extra_surface = Some((Some(surface_idx), plane));
                    }

                    if extra_surface.is_none()
                        && let Some(plane) = find_topological_plane_for_edge_topods(self.tbrep(), edge_idx)
                    {
                        extra_surface = Some((None, plane));
                    }

                    if extra_surface.is_none()
                        && let Some(base_plane) = first_plane
                        && count_plane_face_occurrences_for_line_edge_topods(self.tbrep(), edge_idx) >= 2
                    {
                        extra_surface = Some((None, base_plane));
                    }

                    if let Some((surface_idx, extra_plane)) = extra_surface
                        && let Some(extra_curve2d) =
                            synthesize_plane_pcurve_for_edge_topods(self.tbrep(), edge_idx, &extra_plane)
                    {
                        let extra_surface_id = if let Some(idx) = surface_idx {
                            self.get_or_write_surface_for_face(idx)
                        } else if let Some(pc0) = pcurves.first() {
                            self.get_or_write_surface_for_face(pc0.surface_idx)
                        } else {
                            let placement = self.axis2_from_origin_axis(
                                "face_plane_fallback",
                                extra_plane.origin,
                                extra_plane.normal,
                            );
                            self.write_surface(Some(Surface3::Plane(extra_plane)), Some(placement))
                        };
                        let extra_param = self.write_curve2d(Some(extra_curve2d));
                        let extra_def = self.definitional_representation(extra_param);
                        let extra_pc = self.pcurve(extra_surface_id, extra_def);
                        pcurve_ids.push(extra_pc);
                    }
                }

                let use_torus_seam_curve = pcurve_ids.len() >= 2
                    && torus_surface.is_some()
                    && matches!(
                        source_curve,
                        Some(Curve3::Line(_)) | Some(Curve3::Circle(_)) | None
                    );

                let use_periodic_line_seam_curve =
                    pcurve_ids.len() >= 2 && periodic_line_dup_for_seam;

                let use_seam_curve = use_torus_seam_curve || use_periodic_line_seam_curve;

                if use_seam_curve {
                    self.seam_curve(basis_curve_3d, &pcurve_ids)
                } else {
                    self.surface_curve(basis_curve_3d, &pcurve_ids)
                }
            } else {
                // No pcurves from representations — try helpers via topods.
                if let Some((face_tshape_idx, plane)) = find_plane_surface_for_edge_topods(self.tbrep(), edge_idx) {
                    if let Some(curve2d) = synthesize_plane_pcurve_for_edge_topods(self.tbrep(), edge_idx, &plane) {
                            let promote_plane_line = should_promote_plane_line_pcurve_topods(self.tbrep(), edge_idx);
                        let surface_id = self.get_or_write_surface_for_face(face_tshape_idx);
                            let curve2d = if promote_plane_line {
                                plane_line_pcurve_as_bspline(&curve2d)
                            } else {
                                curve2d
                            };
                            let param_curve_id = self.write_curve2d(Some(curve2d));
                        let def_rep = self.definitional_representation(param_curve_id);
                        let pcurve_id = self.pcurve(surface_id, def_rep);
                        let mut pcurve_ids = vec![pcurve_id];

                        if let Some((peer_surface_idx, peer_plane)) =
                            find_peer_plane_surface_for_line_edge_topods(self.tbrep(), edge_idx, &plane)
                            && let Some(peer_curve2d) =
                                synthesize_plane_pcurve_for_edge_topods(self.tbrep(), edge_idx, &peer_plane)
                        {
                                let peer_curve2d = if promote_plane_line {
                                    plane_line_pcurve_as_bspline(&peer_curve2d)
                                } else {
                                    peer_curve2d
                                };
                            let peer_surface_id = if let Some(idx) = peer_surface_idx {
                                self.get_or_write_surface_for_face(idx)
                            } else {
                                surface_id
                            };
                            let peer_param = self.write_curve2d(Some(peer_curve2d));
                            let peer_def = self.definitional_representation(peer_param);
                            let peer_pc = self.pcurve(peer_surface_id, peer_def);
                            pcurve_ids.push(peer_pc);
                        } else if count_plane_face_occurrences_for_line_edge_topods(self.tbrep(), edge_idx) >= 2
                            && let Some(peer_curve2d) =
                                synthesize_plane_pcurve_for_edge_topods(self.tbrep(), edge_idx, &plane)
                        {
                                let peer_curve2d = if promote_plane_line {
                                    plane_line_pcurve_as_bspline(&peer_curve2d)
                                } else {
                                    peer_curve2d
                                };
                            let peer_param = self.write_curve2d(Some(peer_curve2d));
                            let peer_def = self.definitional_representation(peer_param);
                            let peer_pc = self.pcurve(surface_id, peer_def);
                            pcurve_ids.push(peer_pc);
                        }

                        self.surface_curve(basis_curve_3d, &pcurve_ids)
                    } else {
                        basis_curve_3d
                    }
                } else {
                    if let Some(plane) = find_topological_plane_for_edge_topods(self.tbrep(), edge_idx) {
                        if let Some(curve2d) = synthesize_plane_pcurve_for_edge_topods(self.tbrep(), edge_idx, &plane) {
                            let placement = self.axis2_from_origin_axis("face_plane_fallback", plane.origin, plane.normal);
                            let surface_id = self.write_surface(Some(Surface3::Plane(plane)), Some(placement));
                            let param_curve_id = self.write_curve2d(Some(curve2d));
                            let def_rep = self.definitional_representation(param_curve_id);
                            let pcurve_id = self.pcurve(surface_id, def_rep);
                            self.surface_curve(basis_curve_3d, &[pcurve_id])
                        } else {
                            basis_curve_3d
                        }
                    } else {
                        basis_curve_3d
                    }
                }
            }
        };

        self.edge_geometry_ids.insert(edge_idx, curve_id);
        curve_id
    }

    fn write_edge_curve_geometry_by_index(&mut self, brep: &BRep, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_geometry_ids.get(&edge_idx) {
            return *existing;
        }

        let curve_id = if brep.edges.get(edge_idx).is_none() {
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            self.line("edge_line", origin, vec)
        } else {
            let edge = &brep.edges[edge_idx];
            let start_point = brep
                .vertices
                .get(edge.start)
                .map(|v| dvec3_to_array(v.point))
                .unwrap_or([0.0, 0.0, 0.0]);
            let end_point = brep
                .vertices
                .get(edge.end)
                .map(|v| dvec3_to_array(v.point))
                .unwrap_or([0.0, 0.0, 0.0]);
            let mut basis_curve_3d =
                self.write_basis_curve_for_edge(brep, edge_idx, start_point, end_point);

            let source_curve = brep
                .geom
                .edge_curve
                .get(edge_idx)
                .copied()
                .flatten()
                .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

            let pcurves = brep
                .geom
                .edge_pcurves
                .get(edge_idx)
                .cloned()
                .unwrap_or_default();

            let is_single_plane_closed_circle = !self.strict_plane_closed_ellipse_done
                && edge.start == edge.end
                && pcurves.len() == 1
                && matches!(
                    brep.geom.surfaces.get(pcurves[0].surface_idx),
                    Some(Surface3::Plane(_))
                )
                && matches!(source_curve.as_ref(), Some(Curve3::Circle(_)));

            if is_single_plane_closed_circle
                && let Some(Curve3::Circle(circle)) = source_curve.as_ref()
            {
                let placement =
                    self.axis2_from_origin_axis("edge_plane_closed_axis", circle.center, circle.normal);
                let major = circle.radius.max(1e-9);
                let minor = (circle.radius * (1.0 - 1e-9)).max(1e-9);
                basis_curve_3d = self.ellipse("edge_plane_closed_ellipse", placement, major, minor);
                self.strict_plane_closed_ellipse_done = true;
            }

            let torus_surface = if !pcurves.is_empty() {
                let mut torus = None;
                let mut all_torus = true;
                for pc in &pcurves {
                    match brep.geom.surfaces.get(pc.surface_idx) {
                        Some(Surface3::Torus(t)) => {
                            if torus.is_none() {
                                torus = Some(*t);
                            }
                        }
                        _ => {
                            all_torus = false;
                            break;
                        }
                    }
                }
                if all_torus { torus } else { None }
            } else {
                None
            };

            if let Some(torus) = torus_surface {
                let mut major_seam = false;
                if let Some(pc0) = pcurves.first()
                    && let Some(Curve2d::Line(l)) = brep.geom.curve2ds.get(pc0.curve2d_idx).cloned()
                {
                    // On torus pcurves, dominant +u direction corresponds to a major-circle seam.
                    major_seam = l.direction.x.abs() > l.direction.y.abs();
                }

                let start_pt = brep
                    .vertices
                    .get(edge.start)
                    .map(|v| v.point)
                    .unwrap_or(glam::DVec3::ZERO);
                let end_pt = brep
                    .vertices
                    .get(edge.end)
                    .map(|v| v.point)
                    .unwrap_or(glam::DVec3::ZERO);
                let axis = torus.axis.normalize_or_zero();
                if axis.length_squared() >= 1e-18 {
                    if major_seam {
                        let radial_vec = start_pt - torus.center - axis * (start_pt - torus.center).dot(axis);
                        let seam_radius = radial_vec.length().max(1e-9);
                        let placement =
                            self.axis2_from_origin_axis("edge_torus_major_axis", torus.center, axis);
                        basis_curve_3d =
                            self.circle("edge_torus_major_circle", placement, seam_radius);
                    } else {
                        let mid = (start_pt + end_pt) * 0.5;
                        let radial_raw = mid - torus.center - axis * (mid - torus.center).dot(axis);
                        let radial = if radial_raw.length_squared() < 1e-18 {
                            any_perpendicular_dvec3(axis)
                        } else {
                            radial_raw.normalize_or_zero()
                        };
                        let circle_center = torus.center + radial * torus.major_radius;
                        let circle_normal = axis.cross(radial).normalize_or_zero();
                        if circle_normal.length_squared() >= 1e-18 {
                            let placement = self.axis2_from_origin_axis(
                                "edge_torus_seam_axis",
                                circle_center,
                                circle_normal,
                            );
                            basis_curve_3d = self.circle(
                                "edge_torus_seam_circle",
                                placement,
                                torus.minor_radius.max(1e-9),
                            );
                        }
                    }
                }
            }

            if !pcurves.is_empty() {
                let mut pcurve_ids = Vec::new();
                let mut periodic_extra_curve2d: Option<Curve2d> = None;
                let mut periodic_line_dup_for_seam = false;
                if pcurves.len() == 1
                    && let Some(pc0) = pcurves.first() {
                        let is_periodic_u_surface = matches!(
                            brep.geom.surfaces.get(pc0.surface_idx),
                            Some(Surface3::Torus(_))
                                | Some(Surface3::Cylinder(_))
                                | Some(Surface3::Cone(_))
                        );
                        if is_periodic_u_surface
                            && let Some(Curve2d::Line(l0)) = brep.geom.curve2ds.get(pc0.curve2d_idx).cloned() {
                                let eps = 1e-9;
                                let mut l1 = l0;
                                let mut duplicated = false;
                                if l0.direction.x.abs() <= eps && l0.direction.y.abs() > eps {
                                    l1.origin.x = l0.origin.x + 2.0 * std::f64::consts::PI;
                                    duplicated = true;
                                }
                                if duplicated {
                                    periodic_extra_curve2d = Some(Curve2d::Line(l1));
                                    periodic_line_dup_for_seam = true;
                                }
                            }
                    }
                let mut first_plane: Option<rcad_kernel::geom::Plane> = None;
                let mut first_is_line = false;
                for (pc_i, pc) in pcurves.iter().enumerate() {
                    let surface_id = self.get_or_write_surface_id(brep, pc.surface_idx);
                    let mut curve2d = brep.geom.curve2ds.get(pc.curve2d_idx).cloned();
                    let is_cylinder = matches!(
                        brep.geom.surfaces.get(pc.surface_idx),
                        Some(Surface3::Cylinder(_))
                    );
                    if is_cylinder
                        && let Some(Curve2d::Line(mut l)) = curve2d
                    {
                        // Canonicalize periodic-u pcurves on cylinders for OCCT-style seam handling.
                        let eps = 1e-9;
                        if l.direction.y.abs() <= eps && l.direction.x.abs() > eps {
                            l.direction.x = l.direction.x.abs();
                            l.direction.y = 0.0;
                            let two_pi = 2.0 * std::f64::consts::PI;
                            let u = l.origin.x.rem_euclid(two_pi);
                            if u <= 1e-6 || (two_pi - u) <= 1e-6 {
                                l.origin.x = 0.0;
                            } else {
                                l.origin.x = u;
                            }
                        }
                        curve2d = Some(Curve2d::Line(l));
                    }

                    if pcurves.len() == 1 {
                        first_plane = brep
                            .geom
                            .surfaces
                            .get(pc.surface_idx)
                            .and_then(|s| match s {
                                Surface3::Plane(p) => Some(*p),
                                _ => None,
                            });
                        first_is_line = matches!(curve2d, Some(Curve2d::Line(_)) | None);
                    }

                    let param_curve_id = self.write_curve2d(curve2d);
                    let def_rep = self.definitional_representation(param_curve_id);
                    let pcurve_id = self.pcurve(surface_id, def_rep);
                    pcurve_ids.push(pcurve_id);

                    if pc_i == 0
                        && let Some(extra_curve2d) = periodic_extra_curve2d.clone()
                    {
                        let extra_param = self.write_curve2d(Some(extra_curve2d));
                        let extra_def = self.definitional_representation(extra_param);
                        let extra_pc = self.pcurve(surface_id, extra_def);
                        pcurve_ids.push(extra_pc);
                    }
                }

                if pcurve_ids.len() == 1
                    && first_is_line
                {
                    let mut extra_surface: Option<(Option<usize>, rcad_kernel::geom::Plane)> = None;
                    if let Some(base_plane) = first_plane
                        && let Some((peer_surface_idx, extra_plane)) =
                            find_peer_plane_surface_for_line_edge(brep, edge_idx, base_plane)
                    {
                        extra_surface = Some((peer_surface_idx, extra_plane));
                    }

                    if extra_surface.is_none()
                        && let Some((surface_idx, plane)) = find_plane_surface_for_edge(brep, edge_idx)
                    {
                        extra_surface = Some((Some(surface_idx), plane));
                    }

                    if extra_surface.is_none()
                        && let Some(plane) = find_topological_plane_for_edge(brep, edge_idx)
                    {
                        extra_surface = Some((None, plane));
                    }

                    if extra_surface.is_none()
                        && let Some(base_plane) = first_plane
                        && count_plane_face_occurrences_for_line_edge(brep, edge_idx) >= 2
                    {
                        // Fallback for duplicated-topology solids where two planar
                        // faces share a geometric edge but not a shared edge index.
                        // Emit a second PCURVE to match OCCT two-PCURVE pattern.
                        extra_surface = Some((None, base_plane));
                    }

                    if let Some((surface_idx, extra_plane)) = extra_surface
                        && let Some(extra_curve2d) =
                            synthesize_plane_pcurve_for_edge(brep, edge_idx, &extra_plane)
                    {
                        let extra_surface_id = if let Some(idx) = surface_idx {
                            self.get_or_write_surface_id(brep, idx)
                        } else if let Some(pc0) = pcurves.first() {
                            self.get_or_write_surface_id(brep, pc0.surface_idx)
                        } else {
                            let placement = self.axis2_from_origin_axis(
                                "face_plane_fallback",
                                extra_plane.origin,
                                extra_plane.normal,
                            );
                            self.write_surface(Some(Surface3::Plane(extra_plane)), Some(placement))
                        };
                        let extra_param = self.write_curve2d(Some(extra_curve2d));
                        let extra_def = self.definitional_representation(extra_param);
                        let extra_pc = self.pcurve(extra_surface_id, extra_def);
                        pcurve_ids.push(extra_pc);
                    }
                }

                let use_torus_seam_curve = pcurve_ids.len() >= 2
                    && torus_surface.is_some()
                    && matches!(
                        source_curve,
                        Some(Curve3::Line(_)) | Some(Curve3::Circle(_)) | None
                    );

                let use_periodic_line_seam_curve =
                    pcurve_ids.len() >= 2 && periodic_line_dup_for_seam;

                let use_seam_curve = use_torus_seam_curve || use_periodic_line_seam_curve;

                if use_seam_curve {
                    self.seam_curve(basis_curve_3d, &pcurve_ids)
                } else {
                    self.surface_curve(basis_curve_3d, &pcurve_ids)
                }
            } else {
                if let Some((surface_idx, plane)) = find_plane_surface_for_edge(brep, edge_idx) {
                    if let Some(curve2d) = synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane) {
                            let promote_plane_line = should_promote_plane_line_pcurve(brep, edge_idx);
                        let surface_id = self.get_or_write_surface_id(brep, surface_idx);
                            let curve2d = if promote_plane_line {
                                plane_line_pcurve_as_bspline(&curve2d)
                            } else {
                                curve2d
                            };
                            let param_curve_id = self.write_curve2d(Some(curve2d));
                        let def_rep = self.definitional_representation(param_curve_id);
                        let pcurve_id = self.pcurve(surface_id, def_rep);
                        let mut pcurve_ids = vec![pcurve_id];

                        if let Some((peer_surface_idx, peer_plane)) =
                            find_peer_plane_surface_for_line_edge(brep, edge_idx, plane)
                            && let Some(peer_curve2d) =
                                synthesize_plane_pcurve_for_edge(brep, edge_idx, &peer_plane)
                        {
                                let peer_curve2d = if promote_plane_line {
                                    plane_line_pcurve_as_bspline(&peer_curve2d)
                                } else {
                                    peer_curve2d
                                };
                            let peer_surface_id = if let Some(idx) = peer_surface_idx {
                                self.get_or_write_surface_id(brep, idx)
                            } else {
                                surface_id
                            };
                            let peer_param = self.write_curve2d(Some(peer_curve2d));
                            let peer_def = self.definitional_representation(peer_param);
                            let peer_pc = self.pcurve(peer_surface_id, peer_def);
                            pcurve_ids.push(peer_pc);
                        } else if count_plane_face_occurrences_for_line_edge(brep, edge_idx) >= 2
                            && let Some(peer_curve2d) =
                                synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane)
                        {
                            // Same fallback as above for strict alignment on line edges.
                                let peer_curve2d = if promote_plane_line {
                                    plane_line_pcurve_as_bspline(&peer_curve2d)
                                } else {
                                    peer_curve2d
                                };
                            let peer_param = self.write_curve2d(Some(peer_curve2d));
                            let peer_def = self.definitional_representation(peer_param);
                            let peer_pc = self.pcurve(surface_id, peer_def);
                            pcurve_ids.push(peer_pc);
                        }

                        self.surface_curve(basis_curve_3d, &pcurve_ids)
                    } else {
                        basis_curve_3d
                    }
                } else {
                    if let Some(plane) = find_topological_plane_for_edge(brep, edge_idx) {
                        if let Some(curve2d) = synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane) {
                            let placement = self.axis2_from_origin_axis("face_plane_fallback", plane.origin, plane.normal);
                            let surface_id = self.write_surface(Some(Surface3::Plane(plane)), Some(placement));
                            let param_curve_id = self.write_curve2d(Some(curve2d));
                            let def_rep = self.definitional_representation(param_curve_id);
                            let pcurve_id = self.pcurve(surface_id, def_rep);
                            self.surface_curve(basis_curve_3d, &[pcurve_id])
                        } else {
                            basis_curve_3d
                        }
                    } else {
                        basis_curve_3d
                    }
                }
            }
        };

        self.edge_geometry_ids.insert(edge_idx, curve_id);
        curve_id
    }

    fn write_standalone_wire_curve_by_index(&mut self, brep: &BRep, edge_idx: usize) -> u64 {
        if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) {
            if let Some(curve) = brep.geom.curves.get(curve_idx) {
                let Some(edge) = brep.edges.get(edge_idx) else {
                    return self.write_edge_curve_geometry_by_index(brep, edge_idx);
                };
                let Some(start_point) = brep.vertices.get(edge.start).map(|v| v.point) else {
                    return self.write_edge_curve_geometry_by_index(brep, edge_idx);
                };
                let Some(end_point) = brep.vertices.get(edge.end).map(|v| v.point) else {
                    return self.write_edge_curve_geometry_by_index(brep, edge_idx);
                };
                let range = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r);
                // OCCT-like handling for standalone circular edges:
                // do not trust persisted curve parameter zero blindly after
                // import/rebuild. Reconstruct trimming from topological edge
                // endpoints so exported TRIMMED_CURVE orientation matches the
                // edge traversal seen by viewers (e.g. FreeCAD via OCCT).
                if let Curve3::Circle(circle) = curve
                    && let Some(curve_id) = self.write_standalone_circle_trimmed_from_edge(
                        circle,
                        start_point,
                        end_point,
                        range,
                    )
                {
                    return curve_id;
                }

                let start = dvec3_to_array(start_point);
                let end = dvec3_to_array(end_point);
                let basis = self.write_curve3_entity(curve, start, end);
                if let Some([t0, t1]) = range {
                    return self.trimmed_curve("wire_trimmed_curve", basis, t0, t1);
                }
                if matches!(curve, Curve3::Line(_)) {
                    return self.trimmed_curve("wire_line_segment", basis, 0.0, 1.0);
                }
                return basis;
            }
        }

        let Some(edge) = brep.edges.get(edge_idx) else {
            return self.write_edge_curve_geometry_by_index(brep, edge_idx);
        };
        let Some(start) = brep.vertices.get(edge.start).map(|v| v.point) else {
            return self.write_edge_curve_geometry_by_index(brep, edge_idx);
        };
        let Some(end) = brep.vertices.get(edge.end).map(|v| v.point) else {
            return self.write_edge_curve_geometry_by_index(brep, edge_idx);
        };

        let delta = end - start;
        let length = delta.length();
        if length <= 1e-12 {
            return self.write_edge_curve_geometry_by_index(brep, edge_idx);
        }
        let origin = self.cartesian_point("wire_origin", dvec3_to_array(start));
        let dir = self.direction("wire_dir", normalize(dvec3_to_array(delta)));
        let vec = self.vector("wire_vec", dir, length);
        let basis = self.line("wire_line", origin, vec);
        self.trimmed_curve("wire_segment", basis, 0.0, 1.0)
    }

    fn write_standalone_wire_curve_by_index_topods(&mut self, edge_idx: usize) -> u64 {
        let tbrep = self.tbrep();
        let Some(ed) = tbrep.tshapes.get(edge_idx).and_then(|ts| {
            if let topods::TShape::Edge(e) = &**ts { Some(e.clone()) } else { None }
        }) else {
            return self.write_edge_curve_geometry_by_index_topods(edge_idx);
        };
        let start_point = tbrep.tshapes.get(ed.first.index).and_then(|ts| {
            if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
        }).unwrap_or_default();
        let end_point = tbrep.tshapes.get(ed.last.index).and_then(|ts| {
            if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
        }).unwrap_or_default();
        let range = ed.range;
        // OCCT-like handling for standalone circular edges:
        // do not trust persisted curve parameter zero blindly after
        // import/rebuild. Reconstruct trimming from topological edge
        // endpoints so exported TRIMMED_CURVE orientation matches the
        // edge traversal seen by viewers (e.g. FreeCAD via OCCT).
        if let Some(Curve3::Circle(circle)) = &ed.curve
            && let Some(curve_id) = self.write_standalone_circle_trimmed_from_edge(
                circle,
                start_point,
                end_point,
                Some(range),
            )
        {
            return curve_id;
        }

        let start = dvec3_to_array(start_point);
        let end = dvec3_to_array(end_point);
        if let Some(curve) = &ed.curve {
            let basis = self.write_curve3_entity(curve, start, end);
            return self.trimmed_curve("wire_trimmed_curve", basis, range[0], range[1]);
        }

        let delta = end_point - start_point;
        let length = delta.length();
        if length <= 1e-12 {
            return self.write_edge_curve_geometry_by_index_topods(edge_idx);
        }
        let origin = self.cartesian_point("wire_origin", dvec3_to_array(start_point));
        let dir = self.direction("wire_dir", normalize(dvec3_to_array(delta)));
        let vec = self.vector("wire_vec", dir, length);
        let basis = self.line("wire_line", origin, vec);
        self.trimmed_curve("wire_segment", basis, 0.0, 1.0)
    }

    fn write_standalone_circle_trimmed_from_edge(
        &mut self,
        circle: &rcad_kernel::geom::Circle3,
        start: glam::DVec3,
        end: glam::DVec3,
        range: Option<[f64; 2]>,
    ) -> Option<u64> {
        // OCCT strategy notes:
        // 1) Build circle placement with X direction aligned to edge start.
        // 2) Compute minor sweep from start->end in that local frame.
        // 3) Use imported trim span as a hint to choose minor vs major arc.
        // 4) Preserve trim parameter unit convention (degree/radian).
        // 5) Emit `.F.` when traversal is reversed (t1 < t0 style range).
        let normal = circle.normal.normalize_or_zero();
        if normal.length_squared() <= 1e-24 {
            return None;
        }
        let start_dir = (start - circle.center).normalize_or_zero();
        let end_dir = (end - circle.center).normalize_or_zero();
        if start_dir.length_squared() <= 1e-24 || end_dir.length_squared() <= 1e-24 {
            return None;
        }
        let placement =
            self.axis2_from_origin_axis_ref("wire_circle_axis", circle.center, normal, start_dir);
        let basis = self.circle("wire_circle_basis", placement, circle.radius.max(1e-9));

        let dot = start_dir.dot(end_dir).clamp(-1.0, 1.0);
        let mut minor_sweep = normal.dot(start_dir.cross(end_dir)).atan2(dot);
        if minor_sweep < 0.0 {
            minor_sweep += std::f64::consts::TAU;
        }
        if minor_sweep < 1e-12 {
            minor_sweep = 0.0;
        }

        let major_sweep: f64 = std::f64::consts::TAU - minor_sweep;
        let mut sweep: f64 = minor_sweep;
        if let Some(span_hint) = normalize_arc_span_hint(range)
            && (major_sweep - span_hint).abs() + 1e-9 < (minor_sweep - span_hint).abs()
        {
            sweep = major_sweep;
        }
        let sweep_param = if range_looks_degrees(range) {
            sweep.to_degrees()
        } else {
            sweep
        };

        let reverse = range.is_some_and(|[t0, t1]| t1 < t0);
        let (t0, t1) = if reverse {
            (sweep_param, 0.0)
        } else {
            (0.0, sweep_param)
        };
        Some(self.trimmed_curve("wire_trimmed_curve", basis, t0, t1))
    }

    /// Topods-native variant: cache by tshape index (XOR with high bit to avoid collision).
    fn get_or_write_surface_id_topods(&mut self, face_tshape_idx: usize) -> u64 {
        let cache_key = face_tshape_idx | 0x8000_0000_0000_0000;
        if let Some(existing) = self.surface_ids.get(&cache_key) {
            return *existing;
        }
        let surface = self.tbrep().tshapes.get(face_tshape_idx).and_then(|ts| {
            if let topods::TShape::Face(fd) = &**ts { fd.surface.clone() } else { None }
        });
        let sid = self.write_surface(surface, None);
        self.surface_ids.insert(cache_key, sid);
        sid
    }

    /// Returns the STEP id for a surface from GeomStore, writing it if not yet done.
    fn get_or_write_surface_id(&mut self, brep: &BRep, surface_idx: usize) -> u64 {
        if let Some(existing) = self.surface_ids.get(&surface_idx) {
            return *existing;
        }

        let surface = brep.geom.surfaces.get(surface_idx).cloned();
        // Build a placeholder placement and write the surface entity directly.
        let sid = self.write_surface(surface, None);
        self.surface_ids.insert(surface_idx, sid);
        sid
    }

    /// Writes a 2D curve entity (for use inside DEFINITIONAL_REPRESENTATION).
    fn write_curve2d(&mut self, curve2d: Option<Curve2d>) -> u64 {
        match curve2d {
            Some(Curve2d::Trimmed(tc)) => {
                // Unwrap Trimmed and write the inner curve 锟?the range
                // restriction is metadata, not a separate STEP entity.
                self.write_curve2d(Some((*tc.curve).clone()))
            }
            Some(Curve2d::Line(l)) => {
                let p = self.cartesian_point_2d("pc_origin", [l.origin.x, l.origin.y]);
                let dir = self.direction_2d("pc_dir", normalize2([l.direction.x, l.direction.y]));
                let mag = (l.direction.length()).max(1e-9);
                let vec = self.vector("pc_vec", dir, mag);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::Circle(c)) => {
                let p = self.cartesian_point_2d("pc_center", [c.center.x, c.center.y]);
                let axis = self.direction_2d("pc_axis", [0.0, 1.0]);
                let placement = self.axis2_placement_2d("pc_placement", p, axis);
                self.circle("pcurve_circle", placement, c.radius.max(1e-9))
            }
            Some(Curve2d::Ellipse(e)) => {
                let p = self.cartesian_point_2d("pc_center", [e.center.x, e.center.y]);
                let ref_dir =
                    self.direction_2d("pc_major_dir", normalize2([e.major_dir.x, e.major_dir.y]));
                let placement = self.axis2_placement_2d("pc_placement", p, ref_dir);
                self.ellipse(
                    "pcurve_ellipse",
                    placement,
                    e.major_radius.max(1e-9),
                    e.minor_radius.max(1e-9),
                )
            }
            Some(Curve2d::BSpline(bs)) => self.write_bspline_curve2d(&bs.clone()),
            Some(Curve2d::CircleInvolute(_)) => {
                // Involute PCurve: no dedicated STEP 2D involute entity writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::ArchimedeanSpiral(_)) => {
                // Spiral PCurve: no dedicated STEP 2D spiral writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::LogarithmicSpiral(_)) => {
                // Spiral PCurve: no dedicated STEP 2D spiral writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::SineWave(_)) => {
                // Sine-wave PCurve: no dedicated STEP 2D sine-wave writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::Bezier(_)) => {
                // Bezier PCurve: fall back to degenerate line (no Bezier 2D STEP writer yet)
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::Parabola(_)) | Some(Curve2d::Hyperbola(_)) => {
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            None => {
                // No 2D curve available: fall back to a degenerate line at origin
                // (valid STEP, carries no geometric info).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
        }
    }

    fn write_basis_curve_for_edge(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        let curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|v| *v)
            .and_then(|curve_idx| brep.geom.curves.get(curve_idx))
            .cloned();

        match curve {
            Some(c) => self.write_curve3_entity(&c, start_point, end_point),
            None => {
                // No curve geometry: approximate as straight line between endpoints
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
        }
    }

    /// Write a `Curve3` value as a STEP curve entity.  Shared helper used by
    /// `write_basis_curve_for_edge` and the `Offset` case.
    fn write_curve3_entity(
        &mut self,
        curve: &Curve3,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        match curve {
            Curve3::Line(line) => {
                let origin = self.cartesian_point("line_origin", dvec3_to_array(line.origin));
                let dir = normalize(dvec3_to_array(line.direction));
                let dir_id = self.direction("line_dir", dir);
                let len = vector_length([
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ])
                .max(1e-9);
                let vec_id = self.vector("line_vec", dir_id, len);
                self.line("edge_line", origin, vec_id)
            }
            Curve3::Circle(circle) => {
                let circle_normal = canonicalize_axis_sign(circle.normal);
                let placement =
                    self.axis2_from_origin_axis("circle_axis", circle.center, circle_normal);
                self.circle("edge_circle", placement, circle.radius.max(1e-9))
            }
            Curve3::Ellipse(ellipse) => {
                let placement = self.axis2_from_origin_axis_ref(
                    "ellipse_axis",
                    ellipse.center,
                    ellipse.normal,
                    ellipse.major_dir,
                );
                self.ellipse(
                    "edge_ellipse",
                    placement,
                    ellipse.major_radius.max(1e-9),
                    ellipse.minor_radius.max(1e-9),
                )
            }
            Curve3::BSpline(bs) => self.write_bspline_curve(&bs.clone(), start_point, end_point),
            Curve3::Hyperbola(h) => {
                let placement =
                    self.axis2_from_origin_axis_ref("hyp_axis", h.center, h.normal, h.major_dir);
                self.hyperbola(
                    "edge_hyperbola",
                    placement,
                    h.semi_major.max(1e-9),
                    h.semi_minor.max(1e-9),
                )
            }
            Curve3::Parabola(p) => {
                let placement =
                    self.axis2_from_origin_axis_ref("par_axis", p.vertex, p.normal, p.axis_dir);
                self.parabola("edge_parabola", placement, p.focal_param.max(1e-9))
            }
            Curve3::CircularHelix(_) => {
                // No dedicated helix STEP writer yet; export a tiny line fallback.
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
            Curve3::SineWave(_) => {
                // No dedicated STEP sine-wave writer yet; export a straight line fallback.
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
            Curve3::Offset(o) => {
                let basis_id = self.write_curve3_entity(&o.basis, start_point, end_point);
                let dir_id = self.direction("offset_dir", dvec3_to_array(o.offset_dir));
                self.offset_curve_3d(basis_id, o.offset_distance, dir_id)
            }
            Curve3::Bezier(_) => {
                // Approximate as straight line
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
        }
    }

    /// Write a B_SPLINE_CURVE_WITH_KNOTS entity for a BSpline curve.
    /// Falls back to a straight line if the curve has no control points.
    fn write_bspline_curve(
        &mut self,
        bs: &rcad_kernel::geom::BSplineCurve3,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        if bs.control_points.is_empty() {
            let p0 = self.cartesian_point("bs_origin", start_point);
            let delta = [
                end_point[0] - start_point[0],
                end_point[1] - start_point[1],
                end_point[2] - start_point[2],
            ];
            let magnitude = vector_length(delta).max(1e-9);
            let direction = self.direction("bs_dir", normalize(delta));
            let vector = self.vector("bs_vec", direction, magnitude);
            return self.line("bs_fallback_line", p0, vector);
        }

        // Write control points
        let cp_ids: Vec<u64> = bs
            .control_points
            .iter()
            .map(|&p| self.cartesian_point("bs_cp", dvec3_to_array(p)))
            .collect();

        // Compress knot vector into (multiplicities, knot_values)
        let (mults, knots) = compress_knot_vector(&bs.knots);

        // Determine knot type
        let knot_type = ".UNSPECIFIED.";

        // All weights 1.0 ??non-rational (UNIFORM_RATIONAL if not)
        let rational = bs.weights.iter().any(|&w| (w - 1.0).abs() > 1e-8);

        self.b_spline_curve_with_knots(
            "bspline_curve",
            bs.degree,
            &cp_ids,
            knot_type,
            &mults,
            &knots,
            rational,
            &bs.weights,
        )
    }

    /// Write a B_SPLINE_CURVE_WITH_KNOTS entity for a 2D BSpline PCurve.
    fn write_bspline_curve2d(&mut self, bs: &BSplineCurve2) -> u64 {
        if bs.control_points.is_empty() {
            let p = self.cartesian_point_2d("bs2_origin", [0.0, 0.0]);
            let dir = self.direction_2d("bs2_dir", [1.0, 0.0]);
            let vec = self.vector("bs2_vec", dir, 1e-9);
            return self.line("bs2_fallback_line", p, vec);
        }

        let cp_ids: Vec<u64> = bs
            .control_points
            .iter()
            .map(|p| self.cartesian_point_2d("bs2_cp", [p.x, p.y]))
            .collect();

        let (mults, knots) = compress_knot_vector(&bs.knots);
        let rational = bs.weights.iter().any(|&w| (w - 1.0).abs() > 1e-8);

        self.b_spline_curve_with_knots(
            "bspline_curve2d",
            bs.degree,
            &cp_ids,
            ".UNSPECIFIED.",
            &mults,
            &knots,
            rational,
            &bs.weights,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn b_spline_curve_with_knots(
        &mut self,
        name: &str,
        degree: usize,
        control_points: &[u64],
        knot_type: &str,
        multiplicities: &[usize],
        knots: &[f64],
        rational: bool,
        weights: &[f64],
    ) -> u64 {
        let cp_list = control_points
            .iter()
            .map(|id| format!("#{}", id))
            .collect::<Vec<_>>()
            .join(",");
        let mult_list = multiplicities
            .iter()
            .map(|m| m.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let knot_list = knots
            .iter()
            .map(|k| format!("{:.9}", k))
            .collect::<Vec<_>>()
            .join(",");

        if rational {
            // B_SPLINE_CURVE_WITH_KNOTS + RATIONAL_B_SPLINE_CURVE complex entity
            let weight_list = weights
                .iter()
                .map(|w| format!("{:.6}", w))
                .collect::<Vec<_>>()
                .join(",");
            self.push(format!(
                "( B_SPLINE_CURVE('{}',{},({}),{},.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({}),({}),{}) RATIONAL_B_SPLINE_CURVE(({}) )",
                name, degree, cp_list, knot_type,
                mult_list, knot_list, knot_type,
                weight_list
            ))
        } else {
            self.push(format!(
                "B_SPLINE_CURVE_WITH_KNOTS('{}',{},({}),{},.F.,.F.,({}),({}),{})",
                name, degree, cp_list, knot_type, mult_list, knot_list, knot_type
            ))
        }
    }
    fn vertex_point_by_tshape_idx(&mut self, tshape_idx: usize) -> u64 {
        let tbrep = self.tbrep();
        if let Some(existing) = self.vertex_point_ids.get(&tshape_idx) {
            return *existing;
        }
        let point = tbrep.tshapes.get(tshape_idx).and_then(|ts| {
            if let topods::TShape::Vertex(vd) = &**ts { Some(vd.point) } else { None }
        }).map(dvec3_to_array).unwrap_or([0.0, 0.0, 0.0]);
        let cartesian = self.cartesian_point("vertex_point", point);
        let vertex = self.vertex_point("vertex", cartesian);
        self.vertex_point_ids.insert(tshape_idx, vertex);
        vertex
    }

    fn vertex_point_by_index(&mut self, brep: &BRep, vertex_idx: usize) -> u64 {
        if let Some(existing) = self.vertex_point_ids.get(&vertex_idx) {
            return *existing;
        }

        let point = brep
            .vertices
            .get(vertex_idx)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let cartesian = self.cartesian_point("vertex_point", point);
        let vertex = self.vertex_point("vertex", cartesian);
        self.vertex_point_ids.insert(vertex_idx, vertex);
        vertex
    }

    fn application_context(&mut self, name: &str) -> u64 {
        self.push(format!("APPLICATION_CONTEXT('{}')", name))
    }

    fn application_protocol_definition(
        &mut self,
        status: &str,
        schema: &str,
        year: i32,
        context: u64,
    ) -> u64 {
        self.push(format!(
            "APPLICATION_PROTOCOL_DEFINITION('{}','{}',{},#{})",
            status, schema, year, context
        ))
    }

    fn product_context(&mut self, name: &str, frame: u64, discipline: &str) -> u64 {
        self.push(format!(
            "PRODUCT_CONTEXT('{}',#{},'{}')",
            name, frame, discipline
        ))
    }

    fn product(&mut self, id: &str, name: &str, description: &str, contexts: &[u64]) -> u64 {
        self.push(format!(
            "PRODUCT('{}','{}','{}',({}))",
            id,
            name,
            description,
            refs(contexts)
        ))
    }

    fn product_definition_formation(&mut self, id: &str, description: &str, product: u64) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_FORMATION('{}','{}',#{})",
            id, description, product
        ))
    }

    fn product_definition_context(&mut self, name: &str, frame: u64, life_cycle: &str) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_CONTEXT('{}',#{},'{}')",
            name, frame, life_cycle
        ))
    }

    fn product_definition(
        &mut self,
        id: &str,
        description: &str,
        formation: u64,
        frame: u64,
    ) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION('{}','{}',#{},#{})",
            id, description, formation, frame
        ))
    }

    fn product_definition_shape(&mut self, name: &str, description: &str, definition: u64) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_SHAPE('{}','{}',#{})",
            name, description, definition
        ))
    }

    fn shape_definition_representation(&mut self, shape: u64, representation: u64) -> u64 {
        self.push(format!(
            "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
            shape, representation
        ))
    }

    fn length_unit_meter(&mut self) -> u64 {
        self.push("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string())
    }

    fn plane_angle_unit_radian(&mut self) -> u64 {
        self.push("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string())
    }

    fn solid_angle_unit_steradian(&mut self) -> u64 {
        self.push("( NAMED_UNIT(*) SOLID_ANGLE_UNIT() SI_UNIT($,.STERADIAN.) )".to_string())
    }

    fn uncertainty_measure_with_unit_value(&mut self, length_unit: u64, value: f64) -> u64 {
        self.push(format!(
            "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({:.12E}),#{},'distance_accuracy_value','confusion accuracy')",
            value,
            length_unit
        ))
    }

    fn geometric_representation_context(
        &mut self,
        dimension: i32,
        uncertainty: u64,
        units: &[u64],
        name: &str,
        description: &str,
    ) -> u64 {
        self.push(format!(
            "( GEOMETRIC_REPRESENTATION_CONTEXT({}) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT(({})) REPRESENTATION_CONTEXT('{}','{}') )",
            dimension,
            uncertainty,
            refs(units),
            name,
            description
        ))
    }

    fn cartesian_point(&mut self, name: &str, coords: [f64; 3]) -> u64 {
        self.push(format!(
            "CARTESIAN_POINT('{}',({:.9},{:.9},{:.9}))",
            name, coords[0], coords[1], coords[2]
        ))
    }

    fn cartesian_point_2d(&mut self, name: &str, coords: [f64; 2]) -> u64 {
        self.push(format!(
            "CARTESIAN_POINT('{}',({:.9},{:.9}))",
            name, coords[0], coords[1]
        ))
    }

    fn direction(&mut self, name: &str, coords: [f64; 3]) -> u64 {
        self.push(format!(
            "DIRECTION('{}',({:.9},{:.9},{:.9}))",
            name, coords[0], coords[1], coords[2]
        ))
    }

    fn direction_2d(&mut self, name: &str, coords: [f64; 2]) -> u64 {
        self.push(format!(
            "DIRECTION('{}',({:.9},{:.9}))",
            name, coords[0], coords[1]
        ))
    }

    fn axis2_placement_2d(&mut self, name: &str, location: u64, ref_dir: u64) -> u64 {
        self.push(format!(
            "AXIS2_PLACEMENT_2D('{}',#{},#{})",
            name, location, ref_dir
        ))
    }

    fn definitional_representation(&mut self, curve2d_id: u64) -> u64 {
        // A DEFINITIONAL_REPRESENTATION wraps a 2D curve in a 2D context
        // for use in PCURVE entities.
        let ctx = self.push("( GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('2D SPACE','') )".to_string());
        self.push(format!(
            "DEFINITIONAL_REPRESENTATION('',(#{}),#{})",
            curve2d_id, ctx
        ))
    }

    fn pcurve(&mut self, surface: u64, def_rep: u64) -> u64 {
        self.push(format!("PCURVE('',#{},#{})", surface, def_rep))
    }

    fn surface_curve(&mut self, curve3d: u64, pcurves: &[u64]) -> u64 {
        // SURFACE_CURVE references the 3D curve and a list of associated pcurves.
        // The preference flag .PCURVE_S1. means the first pcurve is preferred for
        // intersection calculations.
        self.push(format!(
            "SURFACE_CURVE('',#{},({}),.PCURVE_S1.)",
            curve3d,
            refs(pcurves)
        ))
    }

    fn seam_curve(&mut self, curve3d: u64, pcurves: &[u64]) -> u64 {
        self.push(format!(
            "SEAM_CURVE('',#{},({}),.PCURVE_S1.)",
            curve3d,
            refs(pcurves)
        ))
    }

    fn vector(&mut self, name: &str, direction: u64, magnitude: f64) -> u64 {
        self.push(format!(
            "VECTOR('{}',#{},{:.9})",
            name, direction, magnitude
        ))
    }

    fn axis2_placement_3d(&mut self, name: &str, origin: u64, axis: u64, ref_dir: u64) -> u64 {
        self.push(format!(
            "AXIS2_PLACEMENT_3D('{}',#{},#{},#{})",
            name, origin, axis, ref_dir
        ))
    }

    fn line(&mut self, name: &str, origin: u64, vector: u64) -> u64 {
        self.push(format!("LINE('{}',#{},#{})", name, origin, vector))
    }

    fn circle(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!("CIRCLE('{}',#{},{:.9})", name, placement, radius))
    }

    fn trimmed_curve(&mut self, name: &str, basis_curve: u64, t0: f64, t1: f64) -> u64 {
        let (trim_start, trim_end, sense) = if t1 < t0 {
            (t1, t0, ".F.")
        } else {
            (t0, t1, ".T.")
        };
        self.push(format!(
            "TRIMMED_CURVE('{}',#{},(PARAMETER_VALUE({:.9})),(PARAMETER_VALUE({:.9})),{},.PARAMETER.)",
            name, basis_curve, trim_start, trim_end, sense
        ))
    }

    fn ellipse(&mut self, name: &str, placement: u64, major: f64, minor: f64) -> u64 {
        self.push(format!(
            "ELLIPSE('{}',#{},{:.9},{:.9})",
            name, placement, major, minor
        ))
    }

    fn hyperbola(&mut self, name: &str, placement: u64, semi_major: f64, semi_minor: f64) -> u64 {
        self.push(format!(
            "HYPERBOLA('{}',#{},{:.9},{:.9})",
            name, placement, semi_major, semi_minor
        ))
    }

    fn parabola(&mut self, name: &str, placement: u64, focal_param: f64) -> u64 {
        self.push(format!(
            "PARABOLA('{}',#{},{:.9})",
            name, placement, focal_param
        ))
    }

    fn offset_curve_3d(&mut self, basis: u64, dist: f64, dir: u64) -> u64 {
        self.push(format!(
            "OFFSET_CURVE_3D('',#{},{:.9},#{})",
            basis, dist, dir
        ))
    }

    fn plane(&mut self, name: &str, placement: u64) -> u64 {
        self.push(format!("PLANE('{}',#{})", name, placement))
    }

    fn cylindrical_surface(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!(
            "CYLINDRICAL_SURFACE('{}',#{},{:.9})",
            name, placement, radius
        ))
    }

    fn spherical_surface(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!(
            "SPHERICAL_SURFACE('{}',#{},{:.9})",
            name, placement, radius
        ))
    }

    fn conical_surface(
        &mut self,
        name: &str,
        placement: u64,
        radius: f64,
        semi_angle_deg: f64,
    ) -> u64 {
        self.push(format!(
            "CONICAL_SURFACE('{}',#{},{:.9},{:.9})",
            name, placement, radius, semi_angle_deg
        ))
    }

    fn toroidal_surface(
        &mut self,
        name: &str,
        placement: u64,
        major_radius: f64,
        minor_radius: f64,
    ) -> u64 {
        self.push(format!(
            "TOROIDAL_SURFACE('{}',#{},{:.9},{:.9})",
            name, placement, major_radius, minor_radius
        ))
    }

    fn axis1_placement(&mut self, name: &str, origin: u64, axis: u64) -> u64 {
        self.push(format!("AXIS1_PLACEMENT('{}',#{},#{})", name, origin, axis))
    }

    fn surface_of_revolution(&mut self, name: &str, swept_curve: u64, axis_placement: u64) -> u64 {
        self.push(format!(
            "SURFACE_OF_REVOLUTION('{}',#{},#{})",
            name, swept_curve, axis_placement
        ))
    }

    fn surface_of_linear_extrusion(
        &mut self,
        name: &str,
        swept_curve: u64,
        extrusion_axis: u64,
    ) -> u64 {
        self.push(format!(
            "SURFACE_OF_LINEAR_EXTRUSION('{}',#{},#{})",
            name, swept_curve, extrusion_axis
        ))
    }

    fn offset_surface(&mut self, name: &str, basis_surface: u64, offset_distance: f64) -> u64 {
        self.push(format!(
            "OFFSET_SURFACE('{}',#{},{:.9},.F.)",
            name, basis_surface, offset_distance
        ))
    }

    fn vertex_point(&mut self, name: &str, point: u64) -> u64 {
        self.push(format!("VERTEX_POINT('{}',#{})", name, point))
    }

    fn edge_curve(
        &mut self,
        name: &str,
        start: u64,
        end: u64,
        curve: u64,
        same_sense: bool,
    ) -> u64 {
        self.push(format!(
            "EDGE_CURVE('{}',#{},#{},#{},{})",
            name,
            start,
            end,
            curve,
            bool_token(same_sense)
        ))
    }

    fn oriented_edge(&mut self, name: &str, edge_curve: u64, orientation: bool) -> u64 {
        self.push(format!(
            "ORIENTED_EDGE('{}',*,*,#{},{})",
            name,
            edge_curve,
            bool_token(orientation)
        ))
    }

    fn edge_loop(&mut self, name: &str, oriented_edges: &[u64]) -> u64 {
        self.push(format!("EDGE_LOOP('{}',({}))", name, refs(oriented_edges)))
    }

    fn vertex_loop(&mut self, name: &str, vertex_point: u64) -> u64 {
        self.push(format!("VERTEX_LOOP('{}',#{})", name, vertex_point))
    }

    fn face_outer_bound(&mut self, name: &str, edge_loop: u64, orientation: bool) -> u64 {
        self.push(format!(
            "FACE_OUTER_BOUND('{}',#{},{})",
            name,
            edge_loop,
            bool_token(orientation)
        ))
    }

    fn face_bound(&mut self, name: &str, edge_loop: u64, orientation: bool) -> u64 {
        self.push(format!(
            "FACE_BOUND('{}',#{},{})",
            name,
            edge_loop,
            bool_token(orientation)
        ))
    }

    fn advanced_face(
        &mut self,
        name: &str,
        bounds: &[u64],
        surface: u64,
        orientation: bool,
    ) -> u64 {
        self.push(format!(
            "ADVANCED_FACE('{}',({}),#{},{})",
            name,
            refs(bounds),
            surface,
            bool_token(orientation)
        ))
    }

    fn open_shell(&mut self, name: &str, faces: &[u64]) -> u64 {
        self.push(format!("OPEN_SHELL('{}',({}))", name, refs(faces)))
    }

    fn closed_shell(&mut self, name: &str, faces: &[u64]) -> u64 {
        self.push(format!("CLOSED_SHELL('{}',({}))", name, refs(faces)))
    }

    fn shell_based_surface_model(&mut self, name: &str, shells: &[u64]) -> u64 {
        self.push(format!(
            "SHELL_BASED_SURFACE_MODEL('{}',({}))",
            name,
            refs(shells)
        ))
    }

    fn manifold_solid_brep(&mut self, name: &str, outer: u64) -> u64 {
        self.push(format!("MANIFOLD_SOLID_BREP('{}',#{})", name, outer))
    }

    fn brep_with_voids(&mut self, name: &str, outer: u64, voids: &[u64]) -> u64 {
        self.push(format!(
            "BREP_WITH_VOIDS('{}',#{},({}))",
            name,
            outer,
            refs(voids),
        ))
    }

    fn compound(&mut self, name: &str, elements: &[u64]) -> u64 {
        self.push(format!("COMPOUND('{}',({}))", name, refs(elements)))
    }

    fn compsolid(&mut self, name: &str, solids: &[u64]) -> u64 {
        self.push(format!("COMPSOLID('{}',({}))", name, refs(solids)))
    }

    fn geometric_curve_set(&mut self, name: &str, curves: &[u64]) -> u64 {
        self.push(format!(
            "GEOMETRIC_CURVE_SET('{}',({}))",
            name,
            refs(curves)
        ))
    }

    fn shape_representation(&mut self, name: &str, items: &[u64], context: u64) -> u64 {
        self.push(format!(
            "SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn manifold_surface_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn advanced_brep_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn geometrically_bounded_wireframe_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn shape_representation_relationship(
        &mut self,
        name: &str,
        description: &str,
        rep_1: u64,
        rep_2: u64,
    ) -> u64 {
        self.push(format!(
            "SHAPE_REPRESENTATION_RELATIONSHIP('{}','{}',#{},#{})",
            name, description, rep_1, rep_2
        ))
    }

    fn item_defined_transformation(
        &mut self,
        name: &str,
        description: &str,
        transform_item_1: u64,
        transform_item_2: u64,
    ) -> u64 {
        self.push(format!(
            "ITEM_DEFINED_TRANSFORMATION('{}','{}',#{},#{})",
            name, description, transform_item_1, transform_item_2
        ))
    }

    fn representation_relationship_with_transformation(
        &mut self,
        name: &str,
        description: &str,
        rep_1: u64,
        rep_2: u64,
        transformation: u64,
    ) -> u64 {
        self.push(format!(
            "( REPRESENTATION_RELATIONSHIP('{}','{}',#{},#{}) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{}) SHAPE_REPRESENTATION_RELATIONSHIP() )",
            name, description, rep_1, rep_2, transformation
        ))
    }

    fn context_dependent_shape_representation(
        &mut self,
        representation_relation: u64,
        represented_product_relation: u64,
    ) -> u64 {
        self.push(format!(
            "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{},#{})",
            representation_relation, represented_product_relation
        ))
    }

    fn next_assembly_usage_occurrence(
        &mut self,
        id: &str,
        name: &str,
        description: &str,
        relating_product_definition: u64,
        related_product_definition: u64,
    ) -> u64 {
        self.push(format!(
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{}','{}',#{},#{},$)",
            id,
            name,
            description,
            relating_product_definition,
            related_product_definition
        ))
    }

    fn product_related_product_category(
        &mut self,
        name: &str,
        description: Option<&str>,
        products: &[u64],
    ) -> u64 {
        let desc = match description {
            Some(v) => format!("'{}'", v),
            None => "$".to_string(),
        };
        self.push(format!(
            "PRODUCT_RELATED_PRODUCT_CATEGORY('{}',{},({}))",
            name,
            desc,
            refs(products)
        ))
    }

    // ?? Color / presentation entities ?????????????????????????????????

    /// Emit STYLED_ITEM + full presentation chain for a single ADVANCED_FACE.
    ///
    /// STEP presentation chain:
    ///   COLOUR_RGB ??FILL_AREA_STYLE_COLOUR ??FILL_AREA_STYLE
    ///   ??SURFACE_STYLE_FILL_AREA ??SURFACE_STYLE_USAGE
    ///   ??SURFACE_SIDE_STYLE ??PRESENTATION_STYLE_ASSIGNMENT
    ///   ??STYLED_ITEM (references the face)
    fn write_face_color(&mut self, face_id: u64, color: Color) {
        let rgb = self.colour_rgb("face_color", color.r, color.g, color.b);
        let fill_color = self.fill_area_style_colour("", rgb);
        let fill_style = self.fill_area_style("", &[fill_color]);
        let fill_area = self.surface_style_fill_area(fill_style);
        let side_style = self.surface_side_style("", &[fill_area]);
        let style_usage = self.surface_style_usage(".BOTH.", side_style);
        let psa = self.presentation_style_assignment(&[style_usage]);
        self.styled_item("color", &[psa], face_id);
    }

    fn write_curve_color(&mut self, curve_id: u64, color: Color, width: f64) {
        let rgb = self.colour_rgb("curve_color", color.r, color.g, color.b);
        let font = self.draughting_pre_defined_curve_font("continuous");
        let curve_style = self.curve_style("", font, width, rgb);
        let psa = self.presentation_style_assignment(&[curve_style]);
        self.styled_item("curve_color", &[psa], curve_id);
    }

    fn colour_rgb(&mut self, name: &str, r: f64, g: f64, b: f64) -> u64 {
        self.push(format!("COLOUR_RGB('{}',{:.6},{:.6},{:.6})", name, r, g, b))
    }

    fn fill_area_style_colour(&mut self, name: &str, colour: u64) -> u64 {
        self.push(format!("FILL_AREA_STYLE_COLOUR('{}',#{})", name, colour))
    }

    fn draughting_pre_defined_curve_font(&mut self, name: &str) -> u64 {
        self.push(format!(
            "DRAUGHTING_PRE_DEFINED_CURVE_FONT('{}')",
            escape_step_string(name)
        ))
    }

    fn curve_style(&mut self, name: &str, font: u64, width: f64, colour: u64) -> u64 {
        self.push(format!(
            "CURVE_STYLE('{}',#{},POSITIVE_LENGTH_MEASURE({:.6}),#{})",
            escape_step_string(name),
            font,
            width,
            colour
        ))
    }

    fn fill_area_style(&mut self, name: &str, styles: &[u64]) -> u64 {
        self.push(format!("FILL_AREA_STYLE('{}',({}))", name, refs(styles)))
    }

    fn surface_style_fill_area(&mut self, style: u64) -> u64 {
        self.push(format!("SURFACE_STYLE_FILL_AREA(#{})", style))
    }

    fn surface_side_style(&mut self, name: &str, styles: &[u64]) -> u64 {
        self.push(format!("SURFACE_SIDE_STYLE('{}',({}))", name, refs(styles)))
    }

    fn surface_style_usage(&mut self, side: &str, style: u64) -> u64 {
        self.push(format!("SURFACE_STYLE_USAGE({},#{})", side, style))
    }

    fn presentation_style_assignment(&mut self, styles: &[u64]) -> u64 {
        self.push(format!("PRESENTATION_STYLE_ASSIGNMENT(({}))", refs(styles)))
    }

    fn styled_item(&mut self, name: &str, styles: &[u64], item: u64) -> u64 {
        self.push(format!(
            "STYLED_ITEM('{}',({},),#{})",
            name,
            refs(styles),
            item
        ))
    }

    fn general_property(&mut self, name: &str, description: &str) -> u64 {
        self.push(format!(
            "GENERAL_PROPERTY('{}','{}',$)",
            escape_step_string(name),
            escape_step_string(description)
        ))
    }

    fn property_definition(&mut self, name: &str, description: &str, reference: u64) -> u64 {
        self.push(format!(
            "PROPERTY_DEFINITION('{}','{}',#{})",
            escape_step_string(name),
            escape_step_string(description),
            reference
        ))
    }

    fn write_ap242_metadata(&mut self, metadata: &StepAp242Metadata) {
        for pdr in &metadata.property_definition_representations {
            let pd = opt_ref_token(pdr.property_definition_id);
            let rep = opt_ref_token(pdr.representation_id);
            self.push(format!("PROPERTY_DEFINITION_REPRESENTATION({},{})", pd, rep));
        }
        for loc in &metadata.dimensional_locations {
            let name = opt_step_string(loc.name.as_deref());
            let desc = opt_step_string(loc.description.as_deref());
            let from_id = opt_ref_token(loc.from_entity_id);
            let to_id = opt_ref_token(loc.to_entity_id);
            self.push(format!(
                "DIMENSIONAL_LOCATION({},{},{},{})",
                name, desc, from_id, to_id
            ));
        }
        for size in &metadata.dimensional_sizes {
            let name = opt_step_string(size.name.as_deref());
            let desc = opt_step_string(size.description.as_deref());
            let shape = opt_ref_token(size.shape_aspect_id);
            self.push(format!("DIMENSIONAL_SIZE({},{},{})", name, desc, shape));
        }
        for tol in &metadata.geometric_tolerances {
            let name = opt_step_string(tol.name.as_deref());
            let desc = opt_step_string(tol.description.as_deref());
            let val = opt_ref_token(tol.value_entity_id);
            let shape = opt_ref_token(tol.shape_aspect_id);
            self.push(format!(
                "GEOMETRIC_TOLERANCE({},{},{},{})",
                name, desc, val, shape
            ));
        }
        for tol in &metadata.geometric_tolerances_with_datum_references {
            let name = opt_step_string(tol.name.as_deref());
            let desc = opt_step_string(tol.description.as_deref());
            let val = opt_ref_token(tol.value_entity_id);
            let shape = opt_ref_token(tol.shape_aspect_id);
            let datum = opt_ref_token(tol.datum_system_id);
            self.push(format!(
                "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE({},{},{},{},{})",
                name, desc, val, shape, datum
            ));
        }
        for datum in &metadata.datums {
            let name = opt_step_string(datum.name.as_deref());
            let desc = opt_step_string(datum.description.as_deref());
            let shape = opt_ref_token(datum.shape_aspect_id);
            self.push(format!("DATUM({},{},{})", name, desc, shape));
        }
        for system in &metadata.datum_systems {
            let name = opt_step_string(system.name.as_deref());
            let desc = opt_step_string(system.description.as_deref());
            let refs = if system.datum_ids.is_empty() {
                "$".to_string()
            } else {
                format!(
                    "({})",
                    system
                        .datum_ids
                        .iter()
                        .map(|id| format!("#{}", id))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            };
            self.push(format!("DATUM_SYSTEM({},{},{})", name, desc, refs));
        }
        for pair in &metadata.kinematic_pairs {
            let entity = sanitize_step_entity_name(&pair.entity_type);
            let name = opt_step_string(pair.name.as_deref());
            let desc = opt_step_string(pair.description.as_deref());
            let refs = if pair.related_entity_ids.is_empty() {
                "$".to_string()
            } else {
                pair
                    .related_entity_ids
                    .iter()
                    .map(|id| format!("#{}", id))
                    .collect::<Vec<_>>()
                    .join(",")
            };
            self.push(format!("{}({},{},{})", entity, name, desc, refs));
        }
        // Write GDT extended entities
        for tol in &metadata.dimensional_tolerances {
            self.write_dimensional_tolerance(tol);
        }
        for val in &metadata.tolerance_values {
            self.write_tolerance_value(val);
        }
        for tol in &metadata.position_tolerances {
            self.write_position_tolerance(tol);
        }
        for tol in &metadata.orientation_tolerances {
            self.write_orientation_tolerance(tol);
        }
        for tol in &metadata.form_tolerances {
            self.write_form_tolerance(tol);
        }
        for tol in &metadata.runout_tolerances {
            self.write_runout_tolerance(tol);
        }
        for tol in &metadata.profile_tolerances {
            self.write_profile_tolerance(tol);
        }
        for frame in &metadata.datum_reference_frames {
            self.write_datum_reference_frame(frame);
        }
        for target in &metadata.datum_targets {
            self.write_datum_target(target);
        }
        for def in &metadata.tolerance_zone_definitions_enhanced {
            self.write_tolerance_zone_definition_enhanced(def);
        }
        // Write view and camera entities
        for view in &metadata.views {
            self.write_view(view);
        }
        for camera in &metadata.cameras {
            self.write_camera_model_d3(camera);
        }
        for volume in &metadata.view_volumes {
            self.write_view_volume(volume);
        }
        // Write annotation entities
        for note in &metadata.notes {
            self.write_note(note);
        }
        for plane in &metadata.annotation_planes {
            self.write_annotation_plane(plane);
        }
        for occurrence in &metadata.annotation_occurrences {
            self.write_annotation_occurrence(occurrence);
        }
        for curve in &metadata.dimension_curves {
            self.write_dimension_curve(curve);
        }
        for symbol in &metadata.terminator_symbols {
            self.write_terminator_symbol(symbol);
        }
        for callout in &metadata.datum_feature_callouts {
            self.write_datum_feature_callout(callout);
        }
    }

    fn write_dimensional_tolerance(&mut self, tol: &StepDimensionalTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let dim_char = opt_ref_token(tol.dimensional_characteristic_id);
        let upper = tol.upper_tolerance.map(|v| format!("{}", v)).unwrap_or_else(|| "$".to_string());
        let lower = tol.lower_tolerance.map(|v| format!("{}", v)).unwrap_or_else(|| "$".to_string());
        let unit = opt_step_string(tol.unit.as_deref());
        self.push(format!(
            "DIMENSIONAL_TOLERANCE({},{},{},{},{},{})",
            name, desc, dim_char, upper, lower, unit
        ));
    }

    fn write_tolerance_value(&mut self, val: &StepToleranceValue) {
        let name = opt_step_string(val.name.as_deref());
        let value = val.value;
        let unit = opt_step_string(val.unit.as_deref());
        self.push(format!(
            "MEASURE_REPRESENTATION_ITEM({},{},{})",
            name, value, unit
        ));
    }

    fn write_position_tolerance(&mut self, tol: &StepPositionTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let projected = if tol.projected { ".T." } else { ".F." };
        let proj_height = tol.projected_height.map(|h| format!("{}", h)).unwrap_or_else(|| "$".to_string());
        self.push(format!(
            "POSITION_TOLERANCE({},{},{},{},{},{},{})",
            name, desc, val, shape, datum, projected, proj_height
        ));
    }

    fn write_orientation_tolerance(&mut self, tol: &StepOrientationTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let entity_name = match tol.orientation_type {
            OrientationToleranceType::Angularity => "ANGULARITY_TOLERANCE",
            OrientationToleranceType::Perpendicularity => "PERPENDICULARITY_TOLERANCE",
            OrientationToleranceType::Parallelism => "PARALLELISM_TOLERANCE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, val, shape, datum
        ));
    }

    fn write_form_tolerance(&mut self, tol: &StepFormTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let entity_name = match tol.form_type {
            FormToleranceType::Flatness => "FLATNESS_TOLERANCE",
            FormToleranceType::Straightness => "STRAIGHTNESS_TOLERANCE",
            FormToleranceType::Circularity => "CIRCULARITY_TOLERANCE",
            FormToleranceType::Cylindricity => "CYLINDRICITY_TOLERANCE",
        };
        self.push(format!("{}({},{},{},{})", entity_name, name, desc, val, shape));
    }

    fn write_runout_tolerance(&mut self, tol: &StepRunoutTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let entity_name = match tol.runout_type {
            RunoutToleranceType::CircularRunout => "CIRCULAR_RUNOUT_TOLERANCE",
            RunoutToleranceType::TotalRunout => "TOTAL_RUNOUT_TOLERANCE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, val, shape, datum
        ));
    }

    fn write_profile_tolerance(&mut self, tol: &StepProfileTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let entity_name = match tol.profile_type {
            ProfileToleranceType::ProfileOfALine => "LINE_PROFILE_TOLERANCE",
            ProfileToleranceType::ProfileOfASurface => "SURFACE_PROFILE_TOLERANCE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, val, shape, datum
        ));
    }

    fn write_datum_reference_frame(&mut self, frame: &StepDatumReferenceFrame) {
        let name = opt_step_string(frame.name.as_deref());
        let desc = opt_step_string(frame.description.as_deref());
        let refs = if frame.datum_system_ids.is_empty() {
            "$".to_string()
        } else {
            format!(
                "({})",
                frame
                    .datum_system_ids
                    .iter()
                    .map(|id| format!("#{}", id))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        self.push(format!("DATUM_REFERENCE_FRAME({},{},{})", name, desc, refs));
    }

    fn write_datum_target(&mut self, target: &StepDatumTarget) {
        let name = opt_step_string(target.name.as_deref());
        let desc = opt_step_string(target.description.as_deref());
        let target_id = opt_step_string(target.target_identifier.as_deref());
        let datum = opt_ref_token(target.datum_id);
        let shape = opt_ref_token(target.shape_aspect_id);
        let entity_name = match target.target_type {
            DatumTargetType::Point => "DATUM_TARGET_POINT",
            DatumTargetType::Line => "DATUM_TARGET_LINE",
            DatumTargetType::Area => "DATUM_TARGET_AREA",
            DatumTargetType::AreaCircle => "DATUM_TARGET_CIRCLE",
            DatumTargetType::AreaRectangle => "DATUM_TARGET_RECTANGLE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, target_id, datum, shape
        ));
    }

    fn write_tolerance_zone_definition_enhanced(&mut self, def: &StepToleranceZoneDefinitionEnhanced) {
        let name = opt_step_string(def.name.as_deref());
        let desc = opt_step_string(def.description.as_deref());
        let zone = opt_ref_token(def.tolerance_zone_id);
        let shape = opt_ref_token(def.shape_aspect_id);
        let defining_shape = opt_ref_token(def.defining_shape_aspect_id);
        self.push(format!(
            "TOLERANCE_ZONE_DEFINITION({},{},{},{},{})",
            name, desc, zone, shape, defining_shape
        ));
    }

    // ?? View and Camera write functions (AP242) ???????????????????????????????

    fn write_view(&mut self, view: &StepView) {
        let name = opt_step_string(view.name.as_deref());
        let desc = opt_step_string(view.description.as_deref());
        let camera = opt_ref_token(view.camera_model_id);
        let view_type = opt_step_string(view.view_type.as_deref());
        self.push(format!(
            "CAMERA_MODEL_D3({},{},{},{})",
            name, desc, camera, view_type
        ));
    }

    fn write_camera_model_d3(&mut self, camera: &StepCameraModelD3) {
        let name = opt_step_string(camera.name.as_deref());
        let view_ref = opt_ref_token(camera.view_reference_system_id);
        let view_volume = opt_ref_token(camera.view_volume_id);
        let perspective = if camera.perspective { ".T." } else { ".F." };
        self.push(format!(
            "CAMERA_MODEL_D3({},{},{},{})",
            name, view_ref, view_volume, perspective
        ));
    }

    fn write_view_volume(&mut self, volume: &StepViewVolume) {
        let name = opt_step_string(volume.name.as_deref());
        let vol_type = match volume.volume_type {
            ViewVolumeType::Orthographic => ".F.",
            ViewVolumeType::Perspective => ".T.",
            ViewVolumeType::Unknown => "$",
        };
        let view_center = volume.view_center
            .map(|v| format!("({},{},{})", v[0], v[1], v[2]))
            .unwrap_or_else(|| "$".to_string());
        let view_plane_dist = volume.view_plane_distance
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "$".to_string());
        let up_dir = volume.up_direction
            .map(|v| format!("({},{},{})", v[0], v[1], v[2]))
            .unwrap_or_else(|| "$".to_string());
        let width = volume.view_window_width
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "$".to_string());
        let height = volume.view_window_height
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "$".to_string());
        self.push(format!(
            "VIEW_VOLUME({},{},{},{},{},{},{})",
            name, vol_type, view_center, view_plane_dist, up_dir, width, height
        ));
    }

    // ?? Annotation write functions (AP242) ????????????????????????????????????

    fn write_note(&mut self, note: &StepNote) {
        let name = opt_step_string(note.name.as_deref());
        let desc = opt_step_string(note.description.as_deref());
        let text = opt_step_string(note.text.as_deref());
        let plane = opt_ref_token(note.annotation_plane_id);
        let geom = opt_ref_token(note.associated_geometry_id);
        self.push(format!(
            "DESCRIPTIVE_REPRESENTATION_ITEM({},{},{},{},{})",
            name, desc, text, plane, geom
        ));
    }

    fn write_annotation_plane(&mut self, plane: &StepAnnotationPlane) {
        let name = opt_step_string(plane.name.as_deref());
        let plane_id = opt_ref_token(plane.plane_id);
        let occurrence = opt_ref_token(plane.annotation_occurrence_id);
        self.push(format!(
            "ANNOTATION_PLANE({},{},{})",
            name, plane_id, occurrence
        ));
    }

    fn write_annotation_occurrence(&mut self, occurrence: &StepAnnotationOccurrence) {
        let name = opt_step_string(occurrence.name.as_deref());
        let style = opt_ref_token(occurrence.style_id);
        let fill = opt_ref_token(occurrence.fill_area_id);
        let shape = opt_ref_token(occurrence.shape_aspect_id);
        self.push(format!(
            "ANNOTATION_OCCURRENCE({},{},{},{})",
            name, style, fill, shape
        ));
    }

    fn write_dimension_curve(&mut self, curve: &StepDimensionCurve) {
        let name = opt_step_string(curve.name.as_deref());
        let curve_id = opt_ref_token(curve.curve_id);
        let plane = opt_ref_token(curve.annotation_plane_id);
        let terminators = if curve.terminator_ids.is_empty() {
            "$".to_string()
        } else {
            format!(
                "({})",
                curve.terminator_ids
                    .iter()
                    .map(|id| format!("#{}", id))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        self.push(format!(
            "DIMENSION_CURVE({},{},{},{})",
            name, curve_id, plane, terminators
        ));
    }

    fn write_terminator_symbol(&mut self, symbol: &StepTerminatorSymbol) {
        let name = opt_step_string(symbol.name.as_deref());
        let curve = opt_ref_token(symbol.annotated_curve_id);
        let term_type = match symbol.terminator_type {
            TerminatorType::Arrow => "FILLED_ARROW",
            TerminatorType::Dot => "FILLED_DOT",
            TerminatorType::OpenArrow => "OPEN_ARROW",
            TerminatorType::ClosedArrow => "FILLED_ARROW",
            TerminatorType::Origin => "ORIGIN_SYMBOL",
            TerminatorType::Unknown => "DIMENSION_CURVE_TERMINATOR",
        };
        self.push(format!(
            "{}({},{})",
            term_type, name, curve
        ));
    }

    fn write_datum_feature_callout(&mut self, callout: &StepDatumFeatureCallout) {
        let name = opt_step_string(callout.name.as_deref());
        let datum_id = opt_step_string(callout.datum_identifier.as_deref());
        let plane = opt_ref_token(callout.annotation_plane_id);
        self.push(format!(
            "DATUM_FEATURE_CALLOUT({},{},{})",
            name, datum_id, plane
        ));
    }

    fn push(&mut self, body: String) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        self.records.push(format!("#{}={};", id, body));
        id
    }
}


// ── Shared utility functions (used by both mod.rs and helpers.rs) ─────────

pub(super) fn esc_step(s: &str) -> String {
    s.replace('\'', "''")
}

pub(super) fn escape_step_string(s: &str) -> String {
    s.replace('\'', "''")
}

pub(super) fn opt_ref_token(v: Option<u64>) -> String {
    v.map(|id| format!("#{}", id))
        .unwrap_or_else(|| "$".to_string())
}

pub(super) fn opt_step_string(s: Option<&str>) -> String {
    s.map(|v| format!("'{}'", escape_step_string(v)))
        .unwrap_or_else(|| "$".to_string())
}

pub(super) fn sanitize_step_entity_name(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    for c in name.chars() {
        if c.is_alphanumeric() || c == '_' || c == '-' || c == '.' {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() { out.push('_'); }
    out
}

// ── Helpers (included inline for same-module visibility) ───────────────────
pub(crate) mod flat;
mod helpers;
use helpers::*; // make pub(super) helpers accessible

#[cfg(test)]
mod tests;
