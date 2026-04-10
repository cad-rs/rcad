use glam::DVec3;
use rcad_kernel::geom::{Surface3, SurfaceEval};

/// 曲面高质量三角化结果。
#[derive(Debug, Clone)]
pub struct SurfaceMesh {
    /// 三角化产生的顶点列表（世界坐标）。
    pub vertices: Vec<DVec3>,
    /// 三角形索引，每个三角形由3个顶点索引组成。
    pub triangles: Vec<[usize; 3]>,
    /// 每个顶点的法向量。
    pub normals: Vec<DVec3>,
}

/// 曲面三角化参数。
#[derive(Debug, Clone)]
pub struct TessellationParams {
    /// 最大弦差（三角形中点到曲面的最大允许距离）。
    /// 较小的值产生更细的网格，推荐范围 0.001~0.1。
    pub chord_tolerance: f64,
    /// 最大角度误差（弧度）。超过此角度的相邻三角形会被进一步细分。
    pub angle_tolerance: f64,
    /// 最小细分步长（UV 空间），防止无限细分。
    pub min_step: f64,
    /// 最大细分步长（UV 空间）。
    pub max_step: f64,
}

impl Default for TessellationParams {
    fn default() -> Self {
        Self {
            chord_tolerance: 0.01,
            angle_tolerance: 0.1,  // ~5.7 degrees
            min_step: 1e-4,
            max_step: 0.5,
        }
    }
}

/// 对参数曲面进行自适应弦差控制三角化。
///
/// 在 UV 参数域上进行自适应细分：
/// 1. 先以均匀网格覆盖 UV 域
/// 2. 对每个四边形检查弦差（三角形中心到真实曲面的距离）
/// 3. 超过 `params.chord_tolerance` 的四边形递归细分
/// 4. 收集所有叶节点三角形
///
/// # 参数
/// - `surface`：要三角化的曲面
/// - `u_range`：UV 域 U 方向范围 [u_min, u_max]
/// - `v_range`：UV 域 V 方向范围 [v_min, v_max]
/// - `params`：三角化参数
pub fn triangulate_surface(
    surface: &Surface3,
    u_range: [f64; 2],
    v_range: [f64; 2],
    params: &TessellationParams,
) -> SurfaceMesh {
    let mut vertices: Vec<DVec3> = Vec::new();
    let mut normals: Vec<DVec3> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // UV 域初始格数（至少 2x2）
    let initial_steps = 4usize;
    let [u0, u1] = u_range;
    let [v0, v1] = v_range;
    let du = (u1 - u0) / initial_steps as f64;
    let dv = (v1 - v0) / initial_steps as f64;

    // 对每个初始四边形进行自适应细分
    for i in 0..initial_steps {
        for j in 0..initial_steps {
            let ua = u0 + i as f64 * du;
            let ub = ua + du;
            let va = v0 + j as f64 * dv;
            let vb = va + dv;

            subdivide_quad(
                surface,
                [ua, ub],
                [va, vb],
                params,
                0,
                &mut vertices,
                &mut normals,
                &mut triangles,
            );
        }
    }

    SurfaceMesh { vertices, triangles, normals }
}

/// 最大递归深度（防止无限细分）。
const MAX_DEPTH: usize = 8;

/// 递归自适应细分一个 UV 四边形。
fn subdivide_quad(
    surface: &Surface3,
    u_range: [f64; 2],
    v_range: [f64; 2],
    params: &TessellationParams,
    depth: usize,
    vertices: &mut Vec<DVec3>,
    normals: &mut Vec<DVec3>,
    triangles: &mut Vec<[usize; 3]>,
) {
    let [u0, u1] = u_range;
    let [v0, v1] = v_range;

    // 计算四角点
    let p00 = surface.point_at(u0, v0);
    let p10 = surface.point_at(u1, v0);
    let p01 = surface.point_at(u0, v1);
    let p11 = surface.point_at(u1, v1);

    let um = (u0 + u1) * 0.5;
    let vm = (v0 + v1) * 0.5;

    // 检查是否需要继续细分
    let should_subdivide = if depth < MAX_DEPTH {
        let step_u = u1 - u0;
        let step_v = v1 - v0;

        // 检查步长是否还能细分
        if step_u < params.min_step * 2.0 && step_v < params.min_step * 2.0 {
            false
        } else {
            // 检查弦差：计算两个三角形的中心点到曲面的距离
            let chord_exceeded = check_chord_tolerance(surface, p00, p10, p11, p01, um, vm, params.chord_tolerance);

            // 检查角度误差（法向量变化）
            let angle_exceeded = depth < MAX_DEPTH / 2 && check_angle_tolerance(surface, u0, u1, v0, v1, params.angle_tolerance);

            chord_exceeded || angle_exceeded
        }
    } else {
        false
    };

    if should_subdivide {
        // 细分为4个子四边形
        subdivide_quad(surface, [u0, um], [v0, vm], params, depth + 1, vertices, normals, triangles);
        subdivide_quad(surface, [um, u1], [v0, vm], params, depth + 1, vertices, normals, triangles);
        subdivide_quad(surface, [u0, um], [vm, v1], params, depth + 1, vertices, normals, triangles);
        subdivide_quad(surface, [um, u1], [vm, v1], params, depth + 1, vertices, normals, triangles);
    } else {
        // 发射两个三角形
        let n = vertices.len();

        // 计算法向量
        let n00 = surface.normal_at(u0, v0);
        let n10 = surface.normal_at(u1, v0);
        let n01 = surface.normal_at(u0, v1);
        let n11 = surface.normal_at(u1, v1);

        // 检查点是否退化（NaN 或 Inf）
        let valid = [p00, p10, p01, p11].iter().all(|p| p.is_finite());
        if !valid {
            return;
        }

        vertices.extend_from_slice(&[p00, p10, p11, p01]);
        normals.extend_from_slice(&[n00, n10, n11, n01]);

        // 选择对角线方向使三角形更均匀
        let d0 = (p11 - p00).length_squared();
        let d1 = (p10 - p01).length_squared();
        if d0 <= d1 {
            // 沿 p00-p11 对角线
            triangles.push([n, n + 1, n + 2]);
            triangles.push([n, n + 2, n + 3]);
        } else {
            // 沿 p10-p01 对角线
            triangles.push([n, n + 1, n + 3]);
            triangles.push([n + 1, n + 2, n + 3]);
        }
    }
}

/// 检查弦差是否超过容差。
/// 计算两个三角形的中心到曲面的近似距离。
fn check_chord_tolerance(
    surface: &Surface3,
    p00: DVec3, p10: DVec3, p11: DVec3, p01: DVec3,
    um: f64, vm: f64,
    tolerance: f64,
) -> bool {
    // 三角形1 (p00, p10, p11) 的中心
    let c1 = (p00 + p10 + p11) / 3.0;
    // 三角形2 (p00, p11, p01) 的中心
    let c2 = (p00 + p11 + p01) / 3.0;

    // 曲面上对应 UV 中心的实际点
    let surf_mid = surface.point_at(um, vm);

    // 检查曲面中点到线性插值中点的距离
    let interp_mid = (p00 + p10 + p11 + p01) / 4.0;
    let chord_err = (surf_mid - interp_mid).length();

    // 也检查各三角形中心处的弦差
    let t1_u = (c1 - p00).length() / (p11 - p00).length().max(1e-10);
    let _ = t1_u; // UV 坐标近似用中点替代
    let chord1 = (surface.point_at(um, vm) - c1).length();
    let chord2 = (surface.point_at(um, vm) - c2).length();

    chord_err > tolerance || chord1 > tolerance || chord2 > tolerance
}

/// 检查角度误差（法向量变化）是否超过容差。
fn check_angle_tolerance(
    surface: &Surface3,
    u0: f64, u1: f64, v0: f64, v1: f64,
    tolerance: f64,
) -> bool {
    let n00 = surface.normal_at(u0, v0);
    let n11 = surface.normal_at(u1, v1);
    let n10 = surface.normal_at(u1, v0);
    let n01 = surface.normal_at(u0, v1);

    // 检查相邻角点的法向量夹角
    for (a, b) in [(n00, n10), (n00, n01), (n11, n10), (n11, n01)] {
        let la = a.length();
        let lb = b.length();
        if la < 0.5 || lb < 0.5 {
            continue;
        }
        let cos_a = (a.dot(b) / (la * lb)).clamp(-1.0, 1.0);
        let angle = cos_a.acos();
        if angle > tolerance {
            return true;
        }
    }
    false
}

/// Ear-clipping triangulation for a simple polygon in 3D.
/// Projects to 2D using the given normal, then runs ear-clipping.
pub fn triangulate_polygon(vertices: &[DVec3], normal: DVec3) -> Vec<[usize; 3]> {
    let n = vertices.len();
    if n < 3 {
        return vec![];
    }
    if n == 3 {
        return vec![[0, 1, 2]];
    }
    if n == 4 {
        return vec![[0, 1, 2], [0, 2, 3]];
    }

    let (u_axis, v_axis) = local_basis(normal);
    let pts_2d: Vec<[f64; 2]> = vertices
        .iter()
        .map(|p| [p.dot(u_axis), p.dot(v_axis)])
        .collect();

    ear_clip(&pts_2d)
}

fn local_basis(normal: DVec3) -> (DVec3, DVec3) {
    let ref_dir = if normal.x.abs() < 0.9 {
        DVec3::X
    } else {
        DVec3::Y
    };
    let u = normal.cross(ref_dir).normalize();
    let v = normal.cross(u).normalize();
    (u, v)
}

fn ear_clip(pts: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let n = pts.len();
    let mut indices: Vec<usize> = (0..n).collect();
    let mut triangles = Vec::new();

    // Ensure CCW winding
    let area = signed_area_2d(pts, &indices);
    if area < 0.0 {
        indices.reverse();
    }

    let mut remaining = indices;
    while remaining.len() > 3 {
        let len = remaining.len();
        let mut ear_found = false;

        for i in 0..len {
            let prev = if i == 0 { len - 1 } else { i - 1 };
            let next = if i == len - 1 { 0 } else { i + 1 };

            let a = remaining[prev];
            let b = remaining[i];
            let c = remaining[next];

            // Check convexity (left turn)
            if cross_2d(pts[a], pts[b], pts[c]) <= 0.0 {
                continue;
            }

            // Check no other vertex inside this triangle
            let mut contains_other = false;
            for j in 0..len {
                if j == prev || j == i || j == next {
                    continue;
                }
                if point_in_triangle_2d(pts[remaining[j]], pts[a], pts[b], pts[c]) {
                    contains_other = true;
                    break;
                }
            }

            if !contains_other {
                triangles.push([a, b, c]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }

        if !ear_found {
            // Degenerate polygon — emit remaining as fan
            for i in 1..remaining.len() - 1 {
                triangles.push([remaining[0], remaining[i], remaining[i + 1]]);
            }
            break;
        }
    }

    if remaining.len() == 3 {
        triangles.push([remaining[0], remaining[1], remaining[2]]);
    }

    triangles
}

fn signed_area_2d(pts: &[[f64; 2]], indices: &[usize]) -> f64 {
    let n = indices.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let a = pts[indices[i]];
        let b = pts[indices[j]];
        area += a[0] * b[1] - b[0] * a[1];
    }
    area * 0.5
}

fn cross_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn point_in_triangle_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let d1 = cross_2d(a, b, p);
    let d2 = cross_2d(b, c, p);
    let d3 = cross_2d(c, a, p);

    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

    !(has_neg && has_pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangulate_triangle() {
        let verts = vec![DVec3::ZERO, DVec3::X, DVec3::Y];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 1);
    }

    #[test]
    fn triangulate_quad() {
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn triangulate_pentagon() {
        let verts = (0..5)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect::<Vec<_>>();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 3);
    }

    #[test]
    fn empty_polygon_returns_no_triangles() {
        let tris = triangulate_polygon(&[], DVec3::Z);
        assert!(tris.is_empty());
    }

    #[test]
    fn two_vertex_polygon_returns_no_triangles() {
        let verts = vec![DVec3::ZERO, DVec3::X];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert!(tris.is_empty());
    }

    #[test]
    fn triangle_count_is_n_minus_2() {
        // A convex n-gon should always yield n-2 triangles.
        for n in 3..=10 {
            let verts: Vec<DVec3> = (0..n)
                .map(|i| {
                    let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                    DVec3::new(a.cos(), a.sin(), 0.0)
                })
                .collect();
            let tris = triangulate_polygon(&verts, DVec3::Z);
            assert_eq!(
                tris.len(),
                n - 2,
                "expected {n}-gon to yield {} triangles, got {}",
                n - 2,
                tris.len()
            );
        }
    }

    #[test]
    fn all_indices_in_bounds() {
        // Every index in the triangulation must be < number of vertices.
        let verts: Vec<DVec3> = (0..7)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 7.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        for tri in &tris {
            for &idx in tri.iter() {
                assert!(idx < verts.len(), "index {idx} out of bounds for {n} vertices", n = verts.len());
            }
        }
    }

    #[test]
    fn clockwise_quad_still_triangulates() {
        // Reversed vertex order (CW) should be handled by sign-flip logic.
        let verts = vec![
            DVec3::new(0.0, 1.0, 0.0), // top-left first (CW)
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
        ];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 2);
    }
}
