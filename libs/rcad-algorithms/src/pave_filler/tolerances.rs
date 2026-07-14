use super::*;

impl<'a> PaveFiller<'a> {
    #[inline]
    pub(crate) fn tol(&self) -> f64 {
        self.ds.fuzzy_tol
    }

    #[inline]
    pub(crate) fn vv_pair_tol(&self, vi: usize, vj: usize) -> f64 {
        self.ds.vertex_tolerance(vi)
            + self.ds.vertex_tolerance(vj)
            + self.tol()
    }

    #[inline]
    pub(crate) fn ve_tol(&self, vi: usize, ei: usize) -> f64 {
        self.ds.vertex_tolerance(vi)
            + self.ds.edge_tolerance(ei)
            + self.tol()
    }

    #[inline]
    pub(crate) fn ee_tol(&self, e1: usize, e2: usize) -> f64 {
        self.ds.edge_tolerance(e1)
            + self.ds.edge_tolerance(e2)
            + self.tol()
    }

    #[inline]
    pub(crate) fn vf_tol(&self, vi: usize, fi: usize) -> f64 {
        self.ds.vertex_tolerance(vi)
            + self.ds.face_tolerance(fi)
            + self.tol()
    }

    #[inline]
    pub(crate) fn ef_tol(&self, ei: usize, fi: usize) -> f64 {
        self.ds.edge_tolerance(ei)
            + self.ds.face_tolerance(fi)
            + self.tol()
    }

    #[inline]
    pub(crate) fn ff_tol(&self, f1: usize, f2: usize) -> f64 {
        self.tol()
            .max(self.ds.face_tolerance(f1))
            .max(self.ds.face_tolerance(f2))
            .max(self.seam_shift_tol)
    }

    /// OCCT: find FF curve indices by face pair
    pub(crate) fn find_face_face_curve_indices(&self, f1: usize, f2: usize) -> Option<Vec<usize>> {
        for ff in &self.ds.interf_ff {
            if ff.f1 == f1 && ff.f2 == f2 {
                return Some(ff.curves.clone());
            }
        }
        None
    }

    /// OCCT: sample face boundary points
    pub(crate) fn sampled_face_boundary_points(&self, face_idx: usize, samples_per_edge: usize) -> Vec<DVec3> {
        let Some(face) = self.ds.faces.get(face_idx) else { return vec![] };
        let mut pts = Vec::new();
        for &ei in &face.boundary_edges {
            if let Some(edge) = self.ds.edges.get(ei) {
                let [t0, t1] = edge.t_range;
                let n = samples_per_edge.max(1);
                for k in 0..=n {
                    let t = t0 + (t1 - t0) * k as f64 / n as f64;
                    let p = edge.curve.point_at(t);
                    if p.is_finite() {
                        pts.push(p);
                    }
                }
            }
        }
        if pts.is_empty() {
            self.ds.face_boundary_points(face_idx)
        } else {
            pts
        }
    }

    /// OCCT: closest point on boundary samples
    pub(crate) fn closest_point_on_boundary_samples(&self, point: DVec3, samples: &[DVec3]) -> DVec3 {
        samples
            .iter()
            .copied()
            .min_by(|a, b| {
                let da = (*a - point).length_squared();
                let db = (*b - point).length_squared();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(point)
    }

    /// OCCT: snap polyline to face boundaries
    pub(crate) fn snap_polyline_endpoints_to_face_boundaries(
        &self,
        chain: &mut Vec<DVec3>,
        f1: usize,
        f2: usize,
    ) {
        if chain.len() < 2 {
            return;
        }

        let boundary_a = self.sampled_face_boundary_points(f1, 12);
        let boundary_b = self.sampled_face_boundary_points(f2, 12);
        if boundary_a.is_empty() || boundary_b.is_empty() {
            return;
        }

        let snap_start_a = self.closest_point_on_boundary_samples(chain[0], &boundary_a);
        let snap_start_b = self.closest_point_on_boundary_samples(chain[0], &boundary_b);
        let snap_end_a = self.closest_point_on_boundary_samples(chain[chain.len() - 1], &boundary_a);
        let snap_end_b = self.closest_point_on_boundary_samples(chain[chain.len() - 1], &boundary_b);

        let choose_better = |orig: DVec3, p1: DVec3, p2: DVec3| {
            let d1 = (p1 - orig).length_squared();
            let d2 = (p2 - orig).length_squared();
            if d1 <= d2 { p1 } else { p2 }
        };

        let start = choose_better(chain[0], snap_start_a, snap_start_b);
        let end = choose_better(chain[chain.len() - 1], snap_end_a, snap_end_b);

        // Only snap if it is a local correction rather than a gross relocation.
        let local_scale = chain
            .windows(2)
            .map(|w| (w[1] - w[0]).length())
            .filter(|d| d.is_finite() && *d > 0.0)
            .fold(f64::INFINITY, f64::min)
            .min(1.0);
        let snap_tol = (local_scale * 4.0)
            .max(TOLERANCE_RETRY_LADDER_COARSE)
            .max(self.ff_tol(f1, f2));

        if (start - chain[0]).length() <= snap_tol {
            chain[0] = start;
        }
        if (end - chain[chain.len() - 1]).length() <= snap_tol {
            let last = chain.len() - 1;
            chain[last] = end;
        }
    }
}
