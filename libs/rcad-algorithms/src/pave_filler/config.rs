use super::*;

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        let n_faces = ds.faces.len();
        let context = IntToolsContext::new(n_faces, TOLERANCE_ABS * 100.0);
        Self {
            ds,
            brep: None,
            face_refs: Vec::new(),
            ic_edge_map: Vec::new(),
            bvh_a: None,
            bvh_b: None,
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
            seam_shift_tol: 0.0,
            // OCCT-aligned: RunParallel (default false)
            run_parallel: false,
            // OCCT-aligned: NonDestructive (default false)
            non_destructive: false,
            // OCCT-aligned: UseOBB (default false)
            use_obb: false,
            // OCCT-aligned: IntTools_Context with FClass2d cache
            context,
            my_arguments: Vec::new(),
            section_attribute: SectionAttribute::default(),
            is_primary: true,
            avoid_build_pcurve: false,
            fpbdone: std::collections::HashMap::new(),
            verts_to_avoid_extension: std::collections::HashSet::new(),
            a_mv_tol: std::collections::HashMap::new(),
            a_dmv_lv: std::collections::HashMap::new(),
            distances: std::collections::HashMap::new(),
        }
    }

    /// Create PaveFiller with BVH acceleration and optional BRep output.
    pub fn with_bvh_and_brep(ds: &'a mut DS, bvh_a: &'a Bvh, bvh_b: &'a Bvh, brep: &'a mut rcad_kernel::topods::BRep) -> Self {
        let total_faces = ds.faces.len();
        let use_bvh = total_faces >= BVH_THRESHOLD;
        let context = IntToolsContext::new(total_faces, TOLERANCE_ABS * 100.0);
        Self {
            ds,
            brep: Some(brep),
            face_refs: Vec::new(),
            ic_edge_map: Vec::new(),
            bvh_a: if use_bvh { Some(bvh_a) } else { None },
            bvh_b: if use_bvh { Some(bvh_b) } else { None },
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
            seam_shift_tol: 0.0,
            run_parallel: false,
            non_destructive: false,
            use_obb: false,
            context,
            my_arguments: Vec::new(),
            section_attribute: SectionAttribute::default(),
            is_primary: true,
            avoid_build_pcurve: false,
            fpbdone: std::collections::HashMap::new(),
            verts_to_avoid_extension: std::collections::HashSet::new(),
            a_mv_tol: std::collections::HashMap::new(),
            a_dmv_lv: std::collections::HashMap::new(),
            distances: std::collections::HashMap::new(),
        }
    }

    /// Create PaveFiller with BVH acceleration.
    pub fn with_bvh(ds: &'a mut DS, bvh_a: &'a Bvh, bvh_b: &'a Bvh) -> Self {
        let total_faces = ds.faces.len();
        let use_bvh = total_faces >= BVH_THRESHOLD;
        let context = IntToolsContext::new(total_faces, TOLERANCE_ABS * 100.0);
        Self {
            ds,
            brep: None,
            face_refs: Vec::new(),
            ic_edge_map: Vec::new(),
            bvh_a: if use_bvh { Some(bvh_a) } else { None },
            bvh_b: if use_bvh { Some(bvh_b) } else { None },
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
            seam_shift_tol: 0.0,
            run_parallel: false,
            non_destructive: false,
            use_obb: false,
            context,
            my_arguments: Vec::new(),
            section_attribute: SectionAttribute::default(),
            is_primary: true,
            avoid_build_pcurve: false,
            fpbdone: std::collections::HashMap::new(),
            verts_to_avoid_extension: std::collections::HashSet::new(),
            a_mv_tol: std::collections::HashMap::new(),
            a_dmv_lv: std::collections::HashMap::new(),
            distances: std::collections::HashMap::new(),
        }
    }

    pub fn configure_glue(&mut self, enable: bool, tolerance: f64) {
        self.use_glue = enable;
        self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
    }

    pub fn configure_glue_adaptive(&mut self, enable: bool, base_tolerance: f64, adaptive: bool) -> f64 {
        if !enable {
            self.use_glue = false;
            return TOLERANCE_ABS;
        }

        self.use_glue = true;

        if !adaptive {
            self.glue_tolerance = base_tolerance.max(TOLERANCE_ABS);
            return self.glue_tolerance;
        }

        // Compute adaptive tolerance based on geometry
        let adaptive_tol = self.compute_adaptive_glue_tolerance(base_tolerance);
        self.glue_tolerance = adaptive_tol;
        adaptive_tol
    }

    pub fn configure_fuzzy(&mut self, fuzzy: f64) {
        self.fuzzy_tolerance = fuzzy.max(0.0);
    }

    pub fn set_run_parallel(&mut self, parallel: bool) {
        self.run_parallel = parallel;
    }

    pub fn set_non_destructive(&mut self, nd: bool) {
        self.non_destructive = nd;
    }

    pub fn set_non_destructive_auto(&mut self) {
        // OCCT: checks if any argument has a locked sub-shape.
        // rcad does not support locked shapes.
        self.non_destructive = false;
    }

    pub fn set_use_obb(&mut self, use_obb: bool) {
        self.use_obb = use_obb;
    }

    /// OCCT-aligned: SetSectionAttribute (BOPAlgo_PaveFiller.hxx L137)
    pub fn set_section_attribute(&mut self, attr: SectionAttribute) {
        self.section_attribute = attr;
    }

    /// OCCT-aligned: SetArguments (BOPAlgo_PaveFiller.hxx L124-127)
    pub fn set_arguments(&mut self, args: Vec<rcad_kernel::BRep>) {
        self.my_arguments = args;
    }

    /// OCCT-aligned: Arguments() const (BOPAlgo_PaveFiller.hxx L133)
    pub fn arguments(&self) -> &[rcad_kernel::BRep] {
        &self.my_arguments
    }

    /// OCCT-aligned: SetAvoidBuildPCurve (BOPAlgo_PaveFiller.hxx L159)
    pub fn set_avoid_build_pcurve(&mut self, flag: bool) {
        self.avoid_build_pcurve = flag;
    }

    /// OCCT-aligned: IsAvoidBuildPCurve() const (BOPAlgo_PaveFiller.hxx L162)
    pub fn is_avoid_build_pcurve(&self) -> bool {
        self.avoid_build_pcurve
    }

    pub fn effective_tolerance(&self, base: f64) -> f64 {
        base.max(self.fuzzy_tolerance)
    }

    fn compute_adaptive_glue_tolerance(&self, base_tolerance: f64) -> f64 {
        let mut min_feature_size = f64::INFINITY;
        let mut min_edge_length = f64::INFINITY;
        let mut min_face_area = f64::INFINITY;

        // Analyze edge lengths
        for edge in &self.ds.edges {
            let p1 = edge.curve.point_at(edge.t_range[0]);
            let p2 = edge.curve.point_at(edge.t_range[1]);
            let length = (p2 - p1).length();
            if length > TOLERANCE_LINEAR_ULTRA_STRICT {
                min_edge_length = min_edge_length.min(length);
            }
        }

        // Analyze face areas (approximate from bounding box)
        for face in &self.ds.faces {
            let pts = self.ds.face_boundary_points(
                self.ds.faces.iter().position(|f| std::ptr::eq(f, face)).unwrap_or(0)
            );
            if pts.len() >= 3 {
                // Compute bounding box diagonal as area proxy
                let mut min_pt = pts[0];
                let mut max_pt = pts[0];
                for p in &pts[1..] {
                    min_pt = min_pt.min(*p);
                    max_pt = max_pt.max(*p);
                }
                let diag = (max_pt - min_pt).length();
                if diag > TOLERANCE_LINEAR_ULTRA_STRICT {
                    min_face_area = min_face_area.min(diag * diag);
                }
            }
        }

        // Use minimum feature size to bound tolerance
        if min_edge_length.is_finite() {
            min_feature_size = min_feature_size.min(min_edge_length);
        }
        if min_face_area.is_finite() {
            min_feature_size = min_feature_size.min(min_face_area.sqrt());
        }

        // Compute adaptive tolerance
        let adaptive_tol = if min_feature_size.is_finite() && min_feature_size > 0.0 {
            // Use a fraction of minimum feature size, but at least base tolerance
            let feature_based = min_feature_size * 0.01;
            base_tolerance.max(feature_based).min(min_feature_size * 0.1)
        } else {
            base_tolerance
        };

        adaptive_tol.max(TOLERANCE_ABS)
    }
}
