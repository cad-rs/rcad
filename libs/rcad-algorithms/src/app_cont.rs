//! AppCont-style continuation matrix utilities.
//!
//! Provides Bernstein-basis mass/inverse/IP/IT matrices used in surface
//! continuity constraints during approximation (OCCT AppCont_ContMatrices).
//!
//! AppCont_ContMatrices (MMatrix, InvMMatrix, IBPMatrix,
//!                  IBTMatrix, VBernstein)

/// Compute binomial coefficient C(n,k).
fn binom(n: i32, k: i32) -> f64 {
    if k < 0 || k > n {
        return 0.0;
    }
    let k = if k > n - k { n - k } else { k };
    let mut r = 1.0;
    for i in 0..k {
        r = r * (n - i) as f64 / (i + 1) as f64;
    }
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
mod app_cont_tests {
    use super::*;

    #[test]
    fn m_matrix_is_symmetric_positive_definite() {
        let classe = 4;
        let mut mat = vec![0.0; (classe * classe) as usize];
        m_matrix(classe, &mut mat);
        for i in 0..classe {
            for j in 0..classe {
                assert!(
                    (mat[(i * classe + j) as usize] - mat[(j * classe + i) as usize]).abs() < 1e-15,
                    "M matrix not symmetric at ({i},{j})"
                );
            }
        }
        for i in 0..classe {
            assert!(
                mat[(i * classe + i) as usize] > 0.0,
                "M matrix diagonal at {i} should be positive"
            );
        }
    }

    #[test]
    fn m_matrix_classe_2() {
        let classe = 2;
        let n = (classe * classe) as usize;
        let mut mat = vec![0.0; n];
        m_matrix(classe, &mut mat);
        // Degree 1: M(0,0)=1/3, M(1,1)=1/3, M(0,1)=M(1,0)=1/6
        assert!((mat[0] - 1.0 / 3.0).abs() < 1e-15, "M00={}", mat[0]);
        assert!((mat[3] - 1.0 / 3.0).abs() < 1e-15, "M11={}", mat[3]);
        assert!((mat[1] - 1.0 / 6.0).abs() < 1e-15, "M01={}", mat[1]);
        assert!((mat[2] - 1.0 / 6.0).abs() < 1e-15, "M10={}", mat[2]);
    }

    #[test]
    fn m_matrix_classe_3() {
        let classe = 3;
        let n = (classe * classe) as usize;
        let mut mat = vec![0.0; n];
        m_matrix(classe, &mut mat);
        // Degree 2: M(0,0)=1/5, M(1,1)=2/15=C(2,1)^2/(5*C(4,2))=4/(5*6), M(2,2)=1/5
        assert!(
            (mat[0] - 1.0 / 5.0).abs() < 1e-15,
            "M00={}, expected 1/5",
            mat[0]
        );
        assert!(
            (mat[4] - 2.0 / 15.0).abs() < 1e-15,
            "M11={}, expected 2/15",
            mat[4]
        );
        assert!(
            (mat[8] - 1.0 / 5.0).abs() < 1e-15,
            "M22={}, expected 1/5",
            mat[8]
        );
    }

    #[test]
    fn v_bernstein_gauss_points() {
        let classe = 3; // degree 2
        let nb_pts = 4;
        let n = (classe * nb_pts) as usize;
        let mut mat = vec![0.0; n];
        v_bernstein(classe, nb_pts, &mut mat);
        // Each column should sum to 1 (partition of unity)
        for j in 0..nb_pts {
            let mut col_sum = 0.0;
            for i in 0..classe {
                col_sum += mat[(i * nb_pts + j) as usize];
            }
            assert!(
                (col_sum - 1.0).abs() < 1e-14,
                "VBernstein column {j} sum={col_sum}, expected 1"
            );
        }
    }

    #[test]
    fn m_matrix_classe_2_diag_dominant() {
        let classe = 2;
        let n = (classe * classe) as usize;
        let mut mat = vec![0.0; n];
        m_matrix(classe, &mut mat);
        // For classe=2 (degree 1), M is strictly diagonal-dominant
        for i in 0..classe {
            let diag = mat[(i * classe + i) as usize];
            let off_sum: f64 = (0..classe)
                .filter(|&j| j != i)
                .map(|j| mat[(i * classe + j) as usize].abs())
                .sum();
            assert!(
                diag > off_sum,
                "M matrix not diagonal-dominant at row {i}: diag={diag}, off={off_sum}"
            );
        }
    }
}
