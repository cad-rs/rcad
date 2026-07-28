//! BRep bounding box computation — delegates to rcad-kernel BndLib.
pub struct BoundingBox;
impl BoundingBox {
    pub fn new() -> rcad_kernel::math::bnd::BndBox { rcad_kernel::math::bnd::BndBox::new() }
}
