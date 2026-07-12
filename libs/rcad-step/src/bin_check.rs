fn main() {
    let mut brep = rcad_kernel::BRep::new();
    let v0 = brep.add_tvertex(glam::DVec3::new(0.0, 0.0, 0.0));
    let v1 = brep.add_tvertex(glam::DVec3::new(1.0, 0.0, 0.0));
    let v2 = brep.add_tvertex(glam::DVec3::new(0.0, 1.0, 0.0));
    let e0 = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
    let e1 = brep.add_tedge(None, v1, v2, [0.0, 1.0]);
    let e2 = brep.add_tedge(None, v2, v0, [0.0, 1.0]);
    let w = brep.add_twire(vec![e0, e1, e2]);
    let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
    let sh = brep.add_tshell(vec![f]);
    brep.add_tsolid(vec![sh]);
    let text = rcad_step::IgesWriter::write_string(&brep);
    println!("IGES output (first 10 lines):");
    for line in text.lines().take(10) {
        println!("  {line}");
    }
    println!("... total lines: {}", text.lines().count());
}
