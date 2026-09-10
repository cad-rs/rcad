//! OCCT AdvApp2Var_SysBase (AdvApp2Var_SysBase.cxx) - the subset consumed by
//! the ApproxF2var / MathBase engines.
//!
//! Encoding notes (architecture difference, annotated at each function):
//! the OCCT primitives operate on raw byte addresses with Fortran workspace
//! offsets (`intptr_t iofset` against a `void* t` reference).  rcad models
//! the workspace as `Vec<f64>` allocations; the offset parameter is kept in
//! the signatures and the C++ semantics (which address is produced) are
//! reproduced, with the Fortran pointer idiom collapsing to index 0 on a
//! freshly allocated buffer.

/// OCCT AdvApp2Var_SysBase::mnfndeb_ (AdvApp2Var_SysBase.cxx L2887-2893) -
/// the Fortran debug level; always 0 in this build (the debug prints of the
/// engines are dead under `ibb >= 2/3`).
pub fn mnfndeb_() -> i32 {
    0
}

/// OCCT AdvApp2Var_SysBase::mgenmsg_ (L2825-2831) - entry debug message;
/// no-op at debug level 0.
pub fn mgenmsg_(_nomprg: &str) -> i32 {
    0
}

/// OCCT AdvApp2Var_SysBase::mgsomsg_ (L2834-2840) - exit debug message;
/// no-op at debug level 0.
pub fn mgsomsg_(_nomprg: &str) -> i32 {
    0
}

/// OCCT AdvApp2Var_SysBase::maermsg_ (L1056-1065) - error message dump;
/// no-op (returns 0, the OCCT body is unconditional).
pub fn maermsg_(_cnompg: &str, _icoder: &mut i32) -> i32 {
    0
}

/// OCCT AdvApp2Var_SysBase::mvriraz_ (L3090-3097) - zero `taille` doubles
/// (the C body memsets `*taille * 8` octets from `adt`).
pub fn mvriraz_(taille: i32, adt: &mut [f64]) {
    let n = (taille * 8) as usize / 8;
    for x in adt.iter_mut().take(n) {
        *x = 0.0;
    }
}

/// OCCT AdvApp2Var_SysBase::miraz_ (L2872-2876) - zero `taille` octets;
/// the consumed callers pass double counts, the byte semantics are kept
/// (taille is octets here, unlike mvriraz_).
pub fn miraz_(taille: i32, adt: &mut [f64]) {
    let n = taille as usize / 8;
    for x in adt.iter_mut().take(n) {
        *x = 0.0;
    }
}

/// OCCT AdvApp2Var_SysBase::mcrfill_ (L2268-2298) - copy `size` OCTETS from
/// `tin` to `tout` (memmove semantics; the OCCT body hand-rolls the
/// overlapping cases around memcpy).  The engines always pass multiples of
/// 8, i.e. whole doubles; the `size / 8` conversion is the rcad encoding of
/// the byte count over f64 storage.
pub fn mcrfill_(size: i32, tin: &[f64], tout: &mut [f64]) {
    let n = size as usize / 8;
    tout.copy_from_slice(&tin[..n]);
}

/// OCCT AdvApp2Var_SysBase::mcrrqst_ (L2504+) - dynamic workspace allocation
/// of `isize` units of `iunit` octets.  Architecture encoding: the C body
/// registers the block and returns `*iofset` such that the caller addresses
/// `t[iofset + i]`; a fresh rcad `Vec<f64>` allocation makes the produced
/// offset 0 and the buffer exactly the requested unit count (the consumed
/// calls are all `iunit = 8` double workspaces).  `iercod` follows the OCCT
/// contract (0 = ok, 3 = allocation refused).
pub fn mcrrqst_(iunit: &mut i32, isize_: &mut i32, t: &mut Vec<f64>, iofset: &mut isize, iercod: &mut i32) {
    if *iunit != 1 && *iunit != 2 && *iunit != 4 && *iunit != 8 {
        *iercod = 2;
        return;
    }
    t.clear();
    t.resize(*isize_ as usize, 0.0);
    *iofset = 0;
    *iercod = 0;
}

/// OCCT AdvApp2Var_SysBase::mcrdelt_ (L1991+) - release of a block acquired
/// by mcrrqst_; the rcad encoding drops the buffer contents (iercod = 0).
pub fn mcrdelt_(iunit: &mut i32, isize_: &mut i32, t: &mut Vec<f64>, iofset: &mut isize, iercod: &mut i32) {
    let _ = (iunit, isize_, iofset);
    t.clear();
    *iercod = 0;
}
