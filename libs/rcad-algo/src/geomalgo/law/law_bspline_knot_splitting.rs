//! OCCT Law_BSplineKnotSplitting (TKGeomAlgo/Law) — 1:1 port of
//! Law_BSplineKnotSplitting.cxx (whole file L26-111) + .hxx members.

use rcad_kernel::math::bspl_lib::max_knot_mult;

use super::law_bspline::LawBSpline;

/// OCCT Law_BSplineKnotSplitting — computes the knot indexes at which the
/// curve must be split to satisfy a given continuity range.
#[derive(Debug, Clone)]
pub struct LawBSplineKnotSplitting {
    split_indexes: Vec<i32>,
}

impl LawBSplineKnotSplitting {
    /// OCCT ctor (Law_BSplineKnotSplitting.cxx L29-89).
    pub fn new(basis_curve: &LawBSpline, continuity_range: i32) -> Self {
        assert!(continuity_range >= 0, "Standard_RangeError");
        let first_index = basis_curve.first_uknot_index();
        let last_index = basis_curve.last_uknot_index();
        let degree = basis_curve.degree();
        let split_indexes = if continuity_range == 0 {
            vec![first_index as i32, last_index as i32]
        } else {
            let mmax = max_knot_mult(
                basis_curve.multiplicities(),
                first_index as i32,
                last_index as i32,
            );
            if (degree as i32 - mmax) >= continuity_range {
                vec![first_index as i32, last_index as i32]
            } else {
                let mults = basis_curve.multiplicities();
                let mut split = vec![0i32; last_index - first_index + 1];
                let mut nb_split = 1usize;
                let mut index = first_index;
                split[nb_split - 1] = index as i32;
                index += 1;
                nb_split += 1;
                while index < last_index {
                    if (degree as i32 - mults[index - 1]) < continuity_range {
                        split[nb_split - 1] = index as i32;
                        nb_split += 1;
                    }
                    index += 1;
                }
                split[nb_split - 1] = index as i32;
                split[..nb_split].to_vec()
            }
        };
        LawBSplineKnotSplitting { split_indexes }
    }

    /// OCCT NbSplits (L91-94).
    pub fn nb_splits(&self) -> usize {
        self.split_indexes.len()
    }

    /// OCCT SplitValue(Index) (L96-102).
    pub fn split_value(&self, index: usize) -> i32 {
        assert!(
            index >= 1 && index <= self.split_indexes.len(),
            "Standard_RangeError"
        );
        self.split_indexes[index - 1]
    }

    /// OCCT Splitting(SplitValues) (L104-110).
    pub fn splitting(&self, split_values: &mut [i32]) {
        for i in 1..=self.split_indexes.len() {
            split_values[i - 1] = self.split_indexes[i - 1];
        }
    }
}
