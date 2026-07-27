use super::*;

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        let n_faces = ds.faces.len();
        let context = IntToolsContext::new(n_faces, TOLERANCE_ABS * 100.0);
        // SAFETY: BOPDS_Iterator borrows ds (lifetime 'a). Transmute to 'static
        // via raw pointer because Rust cannot hold two fields borrowing same data
        // with different mutabilities. PaveFiller<'a> ensures actual lifetime is 'a.
        let my_iterator = {
            let ds_ptr: *const DS = &*ds;
            unsafe {
                std::mem::transmute::<
                    crate::bopds::ds::BOPDS_Iterator<'_>,
                    crate::bopds::ds::BOPDS_Iterator<'static>,
                >(crate::bopds::ds::BOPDS_Iterator::new(unsafe { &*ds_ptr }))
            }
        };
        Self {
            ds,
            my_iterator,
            bvh_a: None,
            bvh_b: None,
            glue: GlueEnum::default(),
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
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
            my_increased_ss: std::collections::HashSet::new(),
            distances: std::collections::HashMap::new(),
            my_report: Report::new(),
            dump_ctx: crate::pipeline_dump::DumpCtx::new_with_module(
                &std::env::var("RCAD_DUMP_GRID").unwrap_or_default(),
                &std::env::var("RCAD_DUMP_CASE").unwrap_or_default(),
                "pf",
            ),
            stop_after: None,
        }
    }

    /// Create PaveFiller with BRep-based BVH (EE/EF/VF pair paths).
    pub fn with_bvh(ds: &'a mut DS, bvh_a: &'a Bvh, bvh_b: &'a Bvh) -> Self {
        let total_faces = ds.faces.len();
        let context = IntToolsContext::new(total_faces, TOLERANCE_ABS * 100.0);
        // SAFETY: see PaveFiller::new()
        let my_iterator = {
            let ds_ptr: *const DS = &*ds;
            unsafe {
                std::mem::transmute::<
                    crate::bopds::ds::BOPDS_Iterator<'_>,
                    crate::bopds::ds::BOPDS_Iterator<'static>,
                >(crate::bopds::ds::BOPDS_Iterator::new(unsafe { &*ds_ptr }))
            }
        };
        Self {
            ds,
            my_iterator,
            bvh_a: if total_faces >= BVH_THRESHOLD {
                Some(bvh_a)
            } else {
                None
            },
            bvh_b: if total_faces >= BVH_THRESHOLD {
                Some(bvh_b)
            } else {
                None
            },
            glue: GlueEnum::default(),
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
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
            my_increased_ss: std::collections::HashSet::new(),
            distances: std::collections::HashMap::new(),
            my_report: Report::new(),
            dump_ctx: crate::pipeline_dump::DumpCtx::new_with_module(
                &std::env::var("RCAD_DUMP_GRID").unwrap_or_default(),
                &std::env::var("RCAD_DUMP_CASE").unwrap_or_default(),
                "pf",
            ),
            stop_after: None,
        }
    }

    // ── continue with existing methods unchanged (configure_glue, configure_fuzzy, etc.) ──

    pub fn configure_glue(&mut self, enable: bool, tolerance: f64) {
        self.glue = if enable {
            GlueEnum::GlueFull
        } else {
            GlueEnum::GlueOff
        };
        self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
    }

    pub fn configure_glue_adaptive(
        &mut self,
        enable: bool,
        base_tolerance: f64,
        adaptive: bool,
    ) -> f64 {
        if !enable {
            self.glue = GlueEnum::GlueOff;
            return TOLERANCE_ABS;
        }
        self.glue = GlueEnum::GlueFull;
        if !adaptive {
            self.glue_tolerance = base_tolerance.max(TOLERANCE_ABS);
            return self.glue_tolerance;
        }
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
        self.non_destructive = false;
    }
    pub fn set_use_obb(&mut self, use_obb: bool) {
        self.use_obb = use_obb;
    }
    pub fn set_section_attribute(&mut self, attr: SectionAttribute) {
        self.section_attribute = attr;
    }
    pub fn set_arguments(&mut self, args: Vec<rcad_kernel::topods::BRep>) {
        self.my_arguments = args;
    }
    pub fn arguments(&self) -> &[rcad_kernel::topods::BRep] {
        &self.my_arguments
    }
    pub fn set_avoid_build_pcurve(&mut self, flag: bool) {
        self.avoid_build_pcurve = flag;
    }
    pub fn is_avoid_build_pcurve(&self) -> bool {
        self.avoid_build_pcurve
    }
    pub fn effective_tolerance(&self, base: f64) -> f64 {
        base.max(self.fuzzy_tolerance)
    }

    fn compute_adaptive_glue_tolerance(&self, base_tolerance: f64) -> f64 {
        let mut min_edge_length = f64::INFINITY;
        let mut min_face_area = f64::INFINITY;
        for edge in &self.ds.edges {
            let p1 = edge.curve.point_at(edge.t_range[0]);
            let p2 = edge.curve.point_at(edge.t_range[1]);
            let length = (p2 - p1).length();
            if length > TOLERANCE_LINEAR_ULTRA_STRICT {
                min_edge_length = min_edge_length.min(length);
            }
        }
        for face in &self.ds.faces {
            let pts = self.ds.face_boundary_points(
                self.ds
                    .faces
                    .iter()
                    .position(|f| std::ptr::eq(f, face))
                    .unwrap_or(0),
            );
            if pts.len() >= 3 {
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
        let mut min_feature_size = f64::INFINITY;
        if min_edge_length.is_finite() {
            min_feature_size = min_feature_size.min(min_edge_length);
        }
        if min_face_area.is_finite() {
            min_feature_size = min_feature_size.min(min_face_area.sqrt());
        }
        if min_feature_size.is_finite() && min_feature_size > 0.0 {
            base_tolerance
                .max(min_feature_size * 0.01)
                .min(min_feature_size * 0.1)
        } else {
            base_tolerance
        }
        .max(TOLERANCE_ABS)
    }
}
