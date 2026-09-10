//! OCCT HLRAlgo static helpers (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo.hxx` (L44-85) + `HLRAlgo.cxx` (L104-292).

use std::sync::LazyLock;

use super::edges_block::MinMaxIndices;

// OCCT HLRAlgo.cxx L106-119: precomputed trig table at PI/14 steps.
static TRIG: LazyLock<[f64; 14]> = LazyLock::new(|| {
    [
        (0. * std::f64::consts::PI / 14.).cos(),
        (0. * std::f64::consts::PI / 14.).sin(),
        (1. * std::f64::consts::PI / 14.).cos(),
        (1. * std::f64::consts::PI / 14.).sin(),
        (2. * std::f64::consts::PI / 14.).cos(),
        (2. * std::f64::consts::PI / 14.).sin(),
        (3. * std::f64::consts::PI / 14.).cos(),
        (3. * std::f64::consts::PI / 14.).sin(),
        (4. * std::f64::consts::PI / 14.).cos(),
        (4. * std::f64::consts::PI / 14.).sin(),
        (5. * std::f64::consts::PI / 14.).cos(),
        (5. * std::f64::consts::PI / 14.).sin(),
        (6. * std::f64::consts::PI / 14.).cos(),
        (6. * std::f64::consts::PI / 14.).sin(),
    ]
});

/// OCCT class HLRAlgo — static MinMax helpers on 16-direction projection
/// boxes.
pub struct HLRAlgo;

impl HLRAlgo {
    /// OCCT UpdateMinMax(x, y, z, Min, Max) — HLRAlgo.cxx L123-158.
    pub fn update_min_max(x: f64, y: f64, z: f64, min: &mut [f64; 16], max: &mut [f64; 16]) {
        let t = &*TRIG;
        let mut d = [0.0f64; 16];
        d[0] = t[0] * x + t[1] * y;
        d[1] = t[1] * x - t[0] * y;
        d[2] = t[2] * x + t[3] * y;
        d[3] = t[3] * x - t[2] * y;
        d[4] = t[4] * x + t[5] * y;
        d[5] = t[5] * x - t[4] * y;
        d[6] = t[6] * x + t[7] * y;
        d[7] = t[7] * x - t[6] * y;
        d[8] = t[8] * x + t[9] * y;
        d[9] = t[9] * x - t[8] * y;
        d[10] = t[10] * x + t[11] * y;
        d[11] = t[11] * x - t[10] * y;
        d[12] = t[12] * x + t[13] * y;
        d[13] = t[13] * x - t[12] * y;
        d[14] = z;
        d[15] = z;

        for i in 0..16 {
            if min[i] > d[i] {
                min[i] = d[i];
            }
            if max[i] < d[i] {
                max[i] = d[i];
            }
        }
    }

    /// OCCT EnlargeMinMax(tol, Min, Max) — HLRAlgo.cxx L162-171.
    pub fn enlarge_min_max(tol: f64, min: &mut [f64; 16], max: &mut [f64; 16]) {
        let mut i = 0;
        while i < 16 {
            min[i] -= tol;
            max[i] += tol;
            i += 1;
        }
    }

    /// OCCT InitMinMax(Big, Min, Max) — HLRAlgo.cxx L175-184.
    pub fn init_min_max(big: f64, min: &mut [f64; 16], max: &mut [f64; 16]) {
        let mut i = 0;
        while i < 16 {
            min[i] = big;
            max[i] = -big;
            i += 1;
        }
    }

    /// OCCT EncodeMinMax(Min, Max, MM) — HLRAlgo.cxx L188-224: two 15-bit
    /// halves packed per 32-bit word, bit-for-bit.
    pub fn encode_min_max(min: &MinMaxIndices, max: &MinMaxIndices, mm: &mut MinMaxIndices) {
        mm.min[0] = min.min[1] & 0x00007fff;
        mm.max[0] = max.min[1] & 0x00007fff;
        mm.min[0] += (min.min[0] & 0x00007fff) << 16;
        mm.max[0] += (max.min[0] & 0x00007fff) << 16;
        mm.min[1] = min.min[3] & 0x00007fff;
        mm.max[1] = max.min[3] & 0x00007fff;
        mm.min[1] += (min.min[2] & 0x00007fff) << 16;
        mm.max[1] += (max.min[2] & 0x00007fff) << 16;
        mm.min[2] = min.min[5] & 0x00007fff;
        mm.max[2] = max.min[5] & 0x00007fff;
        mm.min[2] += (min.min[4] & 0x00007fff) << 16;
        mm.max[2] += (max.min[4] & 0x00007fff) << 16;
        mm.min[3] = min.min[7] & 0x00007fff;
        mm.max[3] = max.min[7] & 0x00007fff;
        mm.min[3] += (min.min[6] & 0x00007fff) << 16;
        mm.max[3] += (max.min[6] & 0x00007fff) << 16;
        mm.min[4] = min.max[1] & 0x00007fff;
        mm.max[4] = max.max[1] & 0x00007fff;
        mm.min[4] += (min.max[0] & 0x00007fff) << 16;
        mm.max[4] += (max.max[0] & 0x00007fff) << 16;
        mm.min[5] = min.max[3] & 0x00007fff;
        mm.max[5] = max.max[3] & 0x00007fff;
        mm.min[5] += (min.max[2] & 0x00007fff) << 16;
        mm.max[5] += (max.max[2] & 0x00007fff) << 16;
        mm.min[6] = min.max[5] & 0x00007fff;
        mm.max[6] = max.max[5] & 0x00007fff;
        mm.min[6] += (min.max[4] & 0x00007fff) << 16;
        mm.max[6] += (max.max[4] & 0x00007fff) << 16;
        mm.min[7] = min.max[7] & 0x00007fff;
        mm.max[7] = max.max[7] & 0x00007fff;
        mm.min[7] += (min.max[6] & 0x00007fff) << 16;
        mm.max[7] += (max.max[6] & 0x00007fff) << 16;
    }

    /// OCCT SizeBox(Min, Max) — HLRAlgo.cxx L228-241.
    pub fn size_box(min: &MinMaxIndices, max: &MinMaxIndices) -> f64 {
        let mut s = (max.min[0] - min.min[0]) as f64;
        for a_i in 1..8 {
            s *= (max.min[a_i] - min.min[a_i]) as f64;
        }
        for a_i in 0..6 {
            s *= (max.max[a_i] - min.max[a_i]) as f64;
        }
        s
    }

    /// OCCT DecodeMinMax(MM, Min, Max) — HLRAlgo.cxx L245-281.
    pub fn decode_min_max(mm: &MinMaxIndices, min: &mut MinMaxIndices, max: &mut MinMaxIndices) {
        min.min[0] = (mm.min[0] & 0x7fff0000) >> 16;
        max.min[0] = (mm.max[0] & 0x7fff0000) >> 16;
        min.min[1] = mm.min[0] & 0x00007fff;
        max.min[1] = mm.max[0] & 0x00007fff;
        min.min[2] = (mm.min[1] & 0x7fff0000) >> 16;
        max.min[2] = (mm.max[1] & 0x7fff0000) >> 16;
        min.min[3] = mm.min[1] & 0x00007fff;
        max.min[3] = mm.max[1] & 0x00007fff;
        min.min[4] = (mm.min[2] & 0x7fff0000) >> 16;
        max.min[4] = (mm.max[2] & 0x7fff0000) >> 16;
        min.min[5] = mm.min[2] & 0x00007fff;
        max.min[5] = mm.max[2] & 0x00007fff;
        min.min[6] = (mm.min[3] & 0x7fff0000) >> 16;
        max.min[6] = (mm.max[3] & 0x7fff0000) >> 16;
        min.min[7] = mm.min[3] & 0x00007fff;
        max.min[7] = mm.max[3] & 0x00007fff;
        min.max[0] = (mm.min[4] & 0x7fff0000) >> 16;
        max.max[0] = (mm.max[4] & 0x7fff0000) >> 16;
        min.max[1] = mm.min[4] & 0x00007fff;
        max.max[1] = mm.max[4] & 0x00007fff;
        min.max[2] = (mm.min[5] & 0x7fff0000) >> 16;
        max.max[2] = (mm.max[5] & 0x7fff0000) >> 16;
        min.max[3] = mm.min[5] & 0x00007fff;
        max.max[3] = mm.max[5] & 0x00007fff;
        min.max[4] = (mm.min[6] & 0x7fff0000) >> 16;
        max.max[4] = (mm.max[6] & 0x7fff0000) >> 16;
        min.max[5] = mm.min[6] & 0x00007fff;
        max.max[5] = mm.max[6] & 0x00007fff;
        min.max[6] = (mm.min[7] & 0x7fff0000) >> 16;
        max.max[6] = (mm.max[7] & 0x7fff0000) >> 16;
        min.max[7] = mm.min[7] & 0x00007fff;
        max.max[7] = mm.max[7] & 0x00007fff;
    }

    /// OCCT CopyMinMax(IMin, IMax, OMin, OMax) — HLRAlgo.hxx L72-79.
    pub fn copy_min_max(
        i_min: &MinMaxIndices,
        i_max: &MinMaxIndices,
        o_min: &mut MinMaxIndices,
        o_max: &mut MinMaxIndices,
    ) {
        *o_min = *i_min;
        *o_max = *i_max;
    }

    /// OCCT AddMinMax(IMin, IMax, OMin, OMax) — HLRAlgo.cxx L285-292.
    pub fn add_min_max(
        i_min: &MinMaxIndices,
        i_max: &MinMaxIndices,
        o_min: &mut MinMaxIndices,
        o_max: &mut MinMaxIndices,
    ) {
        o_min.minimize(i_min);
        o_max.maximize(i_max);
    }
}
