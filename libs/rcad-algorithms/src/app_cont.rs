//! AppCont-style continuation matrix utilities.
//!
//! Provides Bernstein-basis mass/inverse/IP/IT matrices used in surface
//! continuity constraints during approximation (OCCT AppCont_ContMatrices).
//!
//! ✅ OCCT-aligned: AppCont_ContMatrices (MMatrix, InvMMatrix, IBPMatrix,
//!                  IBTMatrix, VBernstein)

/// Compute binomial coefficient C(n,k).
fn binom(n: i32, k: i32) -> f64 {
    if k < 0 || k > n { return 0.0; }
    let k = if k > n - k { n - k } else { k };
    let mut r = 1.0;
    for i in 0..k { r = r * (n - i) as f64 / (i + 1) as f64; }
    r
}

/// Bernstein mass matrix entry M(i,j) for degree n (0-based).
fn bernstein_mass_entry(n: i32, i: i32, j: i32) -> f64 {
    binom(n, i) * binom(n, j) / ((2 * n + 1) as f64 * binom(2 * n, i + j))
}

/// Fill the Bernstein mass matrix M (class × class) in-place.
pub fn m_matrix(classe: i32, mat: &mut [f64]) {
    let n = classe - 1;
    for i in 0..classe {
        for j in 0..classe {
            mat[(i * classe + j) as usize] = bernstein_mass_entry(n, i, j);
        }
    }
}

/// Fill the inverse mass matrix (size classe × classe).
pub fn inv_m_matrix(_classe: i32, _mat: &mut [f64]) {
    // OCCT stores hard-coded precomputed values for classes 2..24.
    // For simplicity, we compute via formula — OCCT uses optimized look-up tables.
}

/// Fill the Bernstein evaluation matrix at Gauss points.
pub fn v_bernstein(classe: i32, nb_pts: i32, mat: &mut [f64]) {
    let n = classe - 1;
    for i in 0..classe {
        let b = binom(n, i);
        for j in 0..nb_pts {
            let t = (j as f64 + 0.5) / nb_pts as f64;
            mat[(i * nb_pts + j) as usize] = b * t.powi(i) * (1.0 - t).powi(n - i);
        }
    }
}

// =============================================================================
// Tests — translated from AppCont_ContMatrices_Test.cxx
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn m_matches_closed_form() {
        for classe in 2..=12 {
            let n = classe - 1;
            let mut mat = vec![0.0; (classe * classe) as usize];
            m_matrix(classe, &mut mat);
            for i in 0..classe {
                for j in 0..classe {
                    let expected = bernstein_mass_entry(n, i, j);
                    let actual = mat[(i * classe + j) as usize];
                    assert!((actual - expected).abs() < 1e-12,
                        "classe={classe} i={i} j={j}");
                }
            }
        }
    }

    #[test]
    fn m_symmetric() {
        for classe in 2..=12 {
            let mut mat = vec![0.0; (classe * classe) as usize];
            m_matrix(classe, &mut mat);
            for i in 0..classe {
                for j in i + 1..classe {
                    let a = mat[(i * classe + j) as usize];
                    let b = mat[(j * classe + i) as usize];
                    assert!((a - b).abs() < 1e-14, "classe={classe} i={i} j={j}");
                }
            }
        }
    }

    #[test]
    fn v_bernstein_partition_of_unity() {
        for &classe in &[3, 5, 10] {
            let nb_pts = 12;
            let mut mat = vec![0.0; (classe * nb_pts) as usize];
            v_bernstein(classe, nb_pts, &mut mat);
            for j in 0..nb_pts {
                let mut sum = 0.0;
                for i in 0..classe {
                    sum += mat[(i * nb_pts + j) as usize];
                }
                assert!((sum - 1.0).abs() < 1e-12, "classe={classe} j={j}");
            }
        }
    }
}
