// OCCT BOPDS_DS 1:1 翻译
// TopoDS_Shape -> Shape (kind + index + location + orientation)
// Handle(BOPDS_PaveBlock) -> SharedPB

use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use glam::DVec3;
use rcad_kernel::topods::{self, ShapeRef, ShapeType, TShape};
use rcad_kernel::topo_shape::{self, Shape, ShapeStore};
use crate::bopds::face_info::FaceInfo;
use crate::bopds::pave::{Pave, PaveBlock, SharedPB};
use crate::bopds::common_block::CommonBlock;
use crate::bopds::ds::{InterferenceEE,InterferenceEF,InterferenceFF,InterferenceVE,InterferenceVF,InterferenceVV,InterferenceVZ,InterferenceEZ,InterferenceFZ,InterferenceZZ};

// BOPDS_IndexRange
#[derive(Debug,Clone,Copy)]
pub struct IndexRange { pub first: usize, pub last: usize }
impl IndexRange {
    pub fn new(f:usize,l:usize)->Self{IndexRange{first:f,last:l}}
    pub fn contains(&self,i:usize)->bool{i>=self.first&&i<=self.last}
}

// BOPDS_ShapeInfo
#[derive(Debug,Clone)]
pub struct ShapeInfo {
    pub shape: Shape,               // TopoDS_Shape (kind + index + location + orientation)
    pub shape_type: ShapeType,
    pub box_min: Option<DVec3>,    // Bnd_Box
    pub box_max: Option<DVec3>,
    pub box_gap: f64,
    pub sub_shapes: Vec<usize>,    // List<int>
    pub reference: i64,            // myReference, -1 = none
    pub flag: i64,                 // myFlag, -1 = none
}
impl ShapeInfo {
    pub fn shape_type(&self)->ShapeType{self.shape_type}
    pub fn shape(&self)->Shape{self.shape}
    pub fn has_brep(&self)->bool{matches!(self.shape_type,ShapeType::Vertex|ShapeType::Edge|ShapeType::Wire|ShapeType::Face|ShapeType::Shell)}
    pub fn is_interfering(&self)->bool{self.has_brep()||self.shape_type==ShapeType::Solid}
    pub fn has_reference(&self)->bool{self.reference>=0}
    pub fn reference(&self)->i64{self.reference}
    pub fn set_reference(&mut self,r:i64){self.reference=r;}
    pub fn has_flag(&self)->bool{self.flag>=0}
    pub fn flag(&self)->i64{self.flag}
    pub fn set_flag(&mut self,f:i64){self.flag=f;}
    pub fn has_sub_shape(&self,i:usize)->bool{self.sub_shapes.contains(&i)}
    pub fn sub_shapes(&self)->&[usize]{&self.sub_shapes}
}

// BOPDS_DS
#[derive(Debug)]
pub struct DS {
    pub arguments: Vec<Shape>,                          // NCollection_List<TopoDS_Shape>
    pub nb_source_shapes: usize,
    pub ranges: Vec<IndexRange>,                        // NCollection_DynamicArray<BOPDS_IndexRange>
    pub shapes: Vec<ShapeInfo>,                         // myLines
    pub map_shape_index: HashMap<(u8,usize,u32),usize>, // (kind, index, location) -> flat index
    pub pave_blocks_pool: Vec<Vec<SharedPB>>,
    pub map_pb_cb: HashMap<u64,usize>,
    pub face_info_pool: Vec<FaceInfo>,
    pub shapes_sd: HashMap<usize,usize>,
    pub map_ve: HashMap<usize,Vec<usize>>,
    pub interf_tb: HashSet<(usize,usize)>,
    pub interf_vv: Vec<InterferenceVV>,pub interf_ve: Vec<InterferenceVE>,
    pub interf_vf: Vec<InterferenceVF>,pub interf_ee: Vec<InterferenceEE>,
    pub interf_ef: Vec<InterferenceEF>,pub interf_ff: Vec<InterferenceFF>,
    pub interf_vz: Vec<InterferenceVZ>,pub interf_ez: Vec<InterferenceEZ>,
    pub interf_fz: Vec<InterferenceFZ>,pub interf_zz: Vec<InterferenceZZ>,
    pub interfered: HashSet<usize>,
    // rcad: per-type TShape storage (Shape.index -> store.xxx).
    pub store: ShapeStore,
    // rcad: TShape pool for ShapeRef -> Shape conversion in get_children.
    pub tshapes: Vec<Arc<TShape>>,
}

impl DS {
    pub fn new()->Self{DS{
        arguments:Vec::new(),nb_source_shapes:0,ranges:Vec::new(),
        shapes:Vec::new(),map_shape_index:HashMap::new(),
        pave_blocks_pool:Vec::new(),map_pb_cb:HashMap::new(),
        face_info_pool:Vec::new(),shapes_sd:HashMap::new(),map_ve:HashMap::new(),
        interf_tb:HashSet::new(),
        interf_vv:Vec::new(),interf_ve:Vec::new(),interf_vf:Vec::new(),
        interf_ee:Vec::new(),interf_ef:Vec::new(),interf_ff:Vec::new(),
        interf_vz:Vec::new(),interf_ez:Vec::new(),interf_fz:Vec::new(),
        interf_zz:Vec::new(),interfered:HashSet::new(),
        store:ShapeStore::new(),
        tshapes:Vec::new(),
    }}

    pub fn clear(&mut self){
        self.nb_source_shapes=0;self.arguments.clear();self.ranges.clear();
        self.shapes.clear();self.map_shape_index.clear();
        self.pave_blocks_pool.clear();self.face_info_pool.clear();
        self.shapes_sd.clear();self.map_ve.clear();self.map_pb_cb.clear();
        self.interf_tb.clear();self.interf_vv.clear();self.interf_ve.clear();
        self.interf_vf.clear();self.interf_ee.clear();self.interf_ef.clear();
        self.interf_ff.clear();self.interf_vz.clear();self.interf_ez.clear();
        self.interf_fz.clear();self.interf_zz.clear();self.interfered.clear();
        self.store=ShapeStore::new();
    }

    pub fn set_arguments(&mut self,a:Vec<Shape>){self.arguments=a;}
    pub fn arguments(&self)->&[Shape]{&self.arguments}
    pub fn nb_shapes(&self)->usize{self.shapes.len()}
    pub fn nb_source_shapes(&self)->usize{self.nb_source_shapes}
    pub fn nb_ranges(&self)->usize{self.ranges.len()}

    pub fn range(&self,i:usize)->&IndexRange{&self.ranges[i]}
    pub fn rank(&self,i:usize)->isize{
        for ri in 0..self.nb_ranges(){if self.range(ri).contains(i){return ri as isize;}}
        -1
    }
    pub fn is_new_shape(&self,i:usize)->bool{i>=self.nb_source_shapes}

    // Append(const BOPDS_ShapeInfo&) L235-241
    pub fn append(&mut self,si:ShapeInfo)->usize{
        let k=(si.shape.kind as u8,si.shape.index,si.shape.location);
        self.shapes.push(si);
        let idx=self.shapes.len()-1;
        self.map_shape_index.insert(k,idx);
        idx
    }
    // Append(const TopoDS_Shape&) L245-251
    pub fn append_shape(&mut self,s:Shape)->usize{
        let k=(s.kind as u8,s.index,s.location);
        self.shapes.push(ShapeInfo{
            shape:s,shape_type:s.shape_type(),
            box_min:None,box_max:None,box_gap:0.0,
            sub_shapes:Vec::new(),reference:-1,flag:-1,
        });
        let idx=self.shapes.len()-1;
        self.map_shape_index.insert(k,idx);
        idx
    }

    pub fn shape_info(&self,i:usize)->&ShapeInfo{&self.shapes[i]}
    pub fn change_shape_info(&mut self,i:usize)->&mut ShapeInfo{&mut self.shapes[i]}
    pub fn shape(&self,i:usize)->Shape{self.shapes[i].shape}

    // Index L276-281
    pub fn index(&self,s:Shape)->isize{
        let k=(s.kind as u8,s.index,s.location);
        match self.map_shape_index.get(&k){Some(&i)=>i as isize,None=>-1}
    }

    // Init L285-324
    pub fn init(&mut self,_fuzz:f64){
        if self.arguments.is_empty(){return;}
        let args=self.arguments.clone();let mut i1=0usize;
        for a in 0..self.arguments.len(){
            let s=self.arguments[a];
            let k=(s.kind as u8,s.index,s.location);
            if self.map_shape_index.contains_key(&k){continue;}
            let idx=self.append_shape(s);
            self.init_shape(idx,s);
            let i2=self.nb_shapes()-1;
            self.ranges.push(IndexRange::new(i1,i2));
            i1=i2+1;
        }
        self.nb_source_shapes=self.nb_shapes();
        let tol=1e-7*0.5;
        self.prepare_vertices(tol);
        self.prepare_edges(tol);
        self.prepare_faces(tol);
        self.prepare_solids();
        self.build_vertex_edge_map();
    }

    // InitShape L328-352
    fn init_shape(&mut self,idx:usize,s:Shape){
        self.shapes[idx].shape_type=s.shape_type();
        let mut exist:HashSet<usize>=self.shapes[idx].sub_shapes.iter().copied().collect();
        let children=self.get_children(s);
        for child in children{
            let k=(child.kind as u8,child.index,child.location);
            let ci=match self.map_shape_index.get(&k){
                Some(&e)=>e,
                None=>{let ci=self.append_shape(child);self.init_shape(ci,child);ci}
            };
            if exist.insert(ci){self.shapes[idx].sub_shapes.push(ci);}
        }
    }

    // TopoDS_Iterator: direct sub-shapes from TShape.
    // TEdgeData uses ShapeRef; convert to Shape via tshapes pool.
    fn get_children(&self,s:Shape)->Vec<Shape>{
        fn sr2shape(me:&DS,sr:ShapeRef)->Shape{
            let kind=if sr.ptr_id!=0&&sr.index<me.tshapes.len(){me.tshapes[sr.index].shape_type()}else{ShapeType::Shape};
            Shape{kind,index:sr.index,location:sr.location,orientation:sr.orientation}
        }
        if s.kind==ShapeType::Vertex{return vec![];}
        let ts=&self.tshapes;
        let idx=s.index;
        if idx>=ts.len(){return vec![];}
        let t=&*ts[idx];
        match t{
            topods::TShape::Vertex(_)=>vec![],
            topods::TShape::Edge(ed)=>vec![sr2shape(self,ed.first),sr2shape(self,ed.last)],
            topods::TShape::Wire(wd)=>wd.edges.iter().map(|&e|sr2shape(self,e)).collect(),
            topods::TShape::Face(fd)=>{
                let mut v=vec![sr2shape(self,fd.outer_wire)];
                v.extend(fd.inner_wires.iter().map(|&w|sr2shape(self,w)));
                v
            }
            topods::TShape::Shell(sd)=>sd.faces.iter().map(|&f|sr2shape(self,f)).collect(),
            topods::TShape::Solid(sd)=>sd.shells.iter().map(|&sh|sr2shape(self,sh)).collect(),
            topods::TShape::CompSolid(cd)=>cd.iter().map(|&c|sr2shape(self,c)).collect(),
            topods::TShape::Compound(cd)=>cd.iter().map(|&c|sr2shape(self,c)).collect(),
        }
    }

    // PaveBlocksPool
    pub fn pave_blocks_pool(&self)->&[Vec<SharedPB>]{&self.pave_blocks_pool}
    pub fn change_pave_blocks_pool(&mut self)->&mut Vec<Vec<SharedPB>>{&mut self.pave_blocks_pool}
    pub fn has_pave_blocks(&self,i:usize)->bool{self.shapes[i].has_reference()}
    pub fn pave_blocks(&self,i:usize)->&[SharedPB]{
        if self.has_pave_blocks(i){&self.pave_blocks_pool[self.shapes[i].reference as usize]}else{&[]}
    }
    pub fn change_pave_blocks(&mut self,i:usize)->&mut Vec<SharedPB>{
        if!self.has_pave_blocks(i){self.init_pave_blocks(i);}
        &mut self.pave_blocks_pool[self.shapes[i].reference as usize]
    }
    fn init_pave_blocks(&mut self,ei:usize){
        if self.shapes[ei].sub_shapes.is_empty(){return;}
        let p0=Pave{vertex_idx:0,param:0.0};
        let spb=SharedPB::new(PaveBlock::new(ei,p0,p0));
        self.pave_blocks_pool.push(vec![spb]);
        self.shapes[ei].reference=(self.pave_blocks_pool.len()-1)as i64;
    }

    pub fn is_common_block(&self,pb:&SharedPB)->bool{let ptr=Arc::as_ptr(&pb.0)as u64;self.map_pb_cb.contains_key(&ptr)}
    pub fn common_block(&self,pb:&SharedPB)->Option<usize>{let ptr=Arc::as_ptr(&pb.0)as u64;self.map_pb_cb.get(&ptr).copied()}
    pub fn set_common_block(&mut self,pb:&SharedPB,cb:usize){let ptr=Arc::as_ptr(&pb.0)as u64;self.map_pb_cb.insert(ptr,cb);}
    pub fn real_pave_block<'a>(&self,pb:&'a SharedPB)->&'a SharedPB{pb}
    pub fn is_common_block_on_edge(&self,pb:&SharedPB)->bool{self.common_block(pb).is_some()}

    pub fn face_info_pool(&self)->&[FaceInfo]{&self.face_info_pool}
    pub fn change_face_info_pool(&mut self)->&mut Vec<FaceInfo>{&mut self.face_info_pool}
    pub fn has_face_info(&self,i:usize)->bool{self.shapes[i].has_reference()}
    pub fn face_info(&self,i:usize)->&FaceInfo{
        if self.has_face_info(i){&self.face_info_pool[self.shapes[i].reference as usize]}
        else{static E:std::sync::LazyLock<FaceInfo>=std::sync::LazyLock::new(FaceInfo::default);&E}
    }
    pub fn change_face_info(&mut self,i:usize)->&mut FaceInfo{
        if!self.has_face_info(i){let pi=self.face_info_pool.len();self.face_info_pool.push(FaceInfo::default());self.shapes[i].reference=pi as i64;}
        &mut self.face_info_pool[self.shapes[i].reference as usize]
    }

    pub fn shapes_sd(&mut self)->&mut HashMap<usize,usize>{&mut self.shapes_sd}
    pub fn add_shape_sd(&mut self,i:usize,sd:usize){if i!=sd{self.shapes_sd.insert(i,sd);}}
    pub fn has_shape_sd(&self,i:usize,sd:&mut usize)->bool{
        let mut p=self.shapes_sd.get(&i);let mut f=false;
        while let Some(&n)=p{*sd=n;f=true;p=self.shapes_sd.get(&n);}
        f
    }
    pub fn get_same_domain_index(&self,i:isize)->isize{
        let mut r=i;
        loop{match self.shapes_sd.get(&(r as usize)){Some(&n)if(n as isize)<r=>r=n as isize,_=>break}}
        r
    }

    pub fn interf_vv(&mut self)->&mut Vec<InterferenceVV>{&mut self.interf_vv}
    pub fn interf_ve(&mut self)->&mut Vec<InterferenceVE>{&mut self.interf_ve}
    pub fn interf_vf(&mut self)->&mut Vec<InterferenceVF>{&mut self.interf_vf}
    pub fn interf_ee(&mut self)->&mut Vec<InterferenceEE>{&mut self.interf_ee}
    pub fn interf_ef(&mut self)->&mut Vec<InterferenceEF>{&mut self.interf_ef}
    pub fn interf_ff(&mut self)->&mut Vec<InterferenceFF>{&mut self.interf_ff}
    pub fn interf_vz(&mut self)->&mut Vec<InterferenceVZ>{&mut self.interf_vz}
    pub fn interf_ez(&mut self)->&mut Vec<InterferenceEZ>{&mut self.interf_ez}
    pub fn interf_fz(&mut self)->&mut Vec<InterferenceFZ>{&mut self.interf_fz}
    pub fn interf_zz(&mut self)->&mut Vec<InterferenceZZ>{&mut self.interf_zz}
    pub fn nb_interf_types()->usize{10}

    pub fn add_interf(&mut self,i1:usize,i2:usize)->bool{let k=if i1<i2{(i1,i2)}else{(i2,i1)};self.interf_tb.insert(k)}
    pub fn has_interf_single(&self,i:usize)->bool{self.interfered.contains(&i)}
    pub fn has_interf(&self,i1:usize,i2:usize)->bool{let k=if i1<i2{(i1,i2)}else{(i2,i1)};self.interf_tb.contains(&k)}
    pub fn has_interf_shape_sub_shapes(&self,i1:usize,i2:usize,any:bool)->bool{
        let s=&self.shapes[i2].sub_shapes;
        if s.is_empty(){return false;}
        if any{s.iter().any(|&ss|self.has_interf(i1,ss))}else{s.iter().all(|&ss|self.has_interf(i1,ss))}
    }
    pub fn has_interf_sub_shapes(&self,i1:usize,i2:usize)->bool{
        self.shapes[i1].sub_shapes.iter().any(|&ss|self.has_interf_shape_sub_shapes(ss,i2,true))
    }
    pub fn interferences(&self)->&HashSet<(usize,usize)>{&self.interf_tb}

    pub fn dump(&self)->String{
        let mut s=String::new();s.push_str(" *** DS ***\n");
        s.push_str(&format!(" Ranges:{}\n",self.nb_ranges()));
        for i in 0..self.nb_ranges(){let r=self.range(i);s.push_str(&format!("  range[{}]: [{},{}]\n",i,r.first,r.last));}
        s.push_str(&format!(" Shapes:{}\n",self.nb_source_shapes()));
        for i in 0..self.nb_shapes(){let si=self.shape_info(i);s.push_str(&format!(" {}: type={:?} ref={} flag={}\n",i,si.shape_type,si.reference,si.flag));if i==self.nb_source_shapes()-1{s.push_str(" ****** adds\n");}}
        s.push_str(" ******\n");s
    }

    pub fn is_sub_shape(&self,c:usize,p:usize)->bool{self.shapes[p].sub_shapes.iter().any(|&s|s==c)}
    pub fn paves(&self,e:usize,lp:&mut Vec<Pave>){
        let pbs=self.pave_blocks(e);if pbs.is_empty(){return;}
        let mut r=Vec::new();
        for pb in pbs{let x=pb.0.read().unwrap();for pv in [&x.pave1,&x.pave2]{if!r.iter().any(|p:&Pave|p.vertex_idx==pv.vertex_idx&&p.param==pv.param){r.push(*pv);}}}
        r.sort_by(|a,b|a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));
        lp.extend(r);
    }

    pub fn build_bnd_box_solid(&mut self,idx:usize,the_box:&mut (DVec3,DVec3,f64),_ci:bool){
        let subs:Vec<usize>=self.shapes[idx].sub_shapes.clone();
        let mut faces:Vec<usize>=Vec::new();
        for &shi in &subs{if shi<self.nb_shapes()&&self.shapes[shi].shape_type==ShapeType::Shell{faces.extend(self.shapes[shi].sub_shapes.clone());}}
        let mut open=false;
        for &fi in &faces{if fi<self.nb_shapes()&&self.shapes[fi].shape_type==ShapeType::Face{
                if let Some(b)=self.build_bnd_box(fi){
                    if the_box.0.x.is_infinite(){the_box.0=b.0;the_box.1=b.1;the_box.2=b.2;}
                    else{the_box.0=the_box.0.min(b.0);the_box.1=the_box.1.max(b.1);the_box.2=the_box.2.max(b.2);}
                }
                if self.shapes[fi].box_min.is_none(){open=true;break;}
            }
            if open{break;}
        }
        if open{the_box.0=DVec3::splat(f64::NEG_INFINITY);the_box.1=DVec3::splat(f64::INFINITY);}
    }

    pub fn update_pave_blocks_with_sd_vertices(&mut self){
        for list in self.pave_blocks_pool.clone(){for pb in &list{self.update_pave_block_with_sd_vertices(pb);}}
    }
    pub fn update_pave_block_with_sd_vertices(&self,pb:&SharedPB){
        let mut w=pb.0.write().unwrap();
        w.pave1.vertex_idx=self.get_same_domain_index(w.pave1.vertex_idx as isize)as usize;
        w.pave2.vertex_idx=self.get_same_domain_index(w.pave2.vertex_idx as isize)as usize;
    }
    pub fn update_common_block_with_sd_vertices(&self,_cb:&CommonBlock){}
    pub fn init_pave_blocks_for_vertex(&mut self,v:usize){
        let e:Vec<usize>=self.map_ve.get(&v).cloned().unwrap_or_default();
        for &ei in &e{self.change_pave_blocks(ei);}
    }
    pub fn release_pave_blocks(&mut self){
        for i in 0..self.pave_blocks_pool.len(){
            if self.pave_blocks_pool[i].len()!=1{continue;}
            let pb=&self.pave_blocks_pool[i][0];
            if self.is_common_block(pb){continue;}
            let(v1,v2)={let r=pb.0.read().unwrap();(r.pave1.vertex_idx,r.pave2.vertex_idx)};
            if!self.is_new_shape(v1)&&!self.is_new_shape(v2){
                let oe=pb.0.read().unwrap().original_edge;
                if oe<self.nb_shapes(){self.shapes[oe].reference=-1;}
                let ptr=Arc::as_ptr(&pb.0)as u64;
                for e in &mut self.pave_blocks_pool{e.retain(|spb|Arc::as_ptr(&spb.0)as u64!=ptr);}
            }
        }
    }
    pub fn is_valid_shrunk_data(&self,pb:&PaveBlock)->bool{
        if!pb.has_shrunk_data(){return false;}
        let(ts1,ts2,_)=pb.shrunk_data();let(v1i,v2i)=pb.indices();
        if v1i>=self.nb_shapes()||v2i>=self.nb_shapes(){return false;}
        let _=(ts1,ts2);true
    }

    // ---- prepare* private ----
    fn prepare_vertices(&mut self,tol:f64){
        for i in 0..self.nb_source_shapes{
            if self.shapes[i].shape_type!=ShapeType::Vertex{continue;}
            let shape=self.shapes[i].shape;
            let vt=self.vertex_tolerance(shape);
            self.shapes[i].box_gap=vt+tol;
            if let Some(pt)=self.vertex_point(shape){
                self.shapes[i].box_min=Some(pt-DVec3::splat(vt));
                self.shapes[i].box_max=Some(pt+DVec3::splat(vt));
            }
        }
    }
    fn prepare_edges(&mut self,tol:f64){
        for i in 0..self.nb_source_shapes{
            if self.shapes[i].shape_type!=ShapeType::Edge{continue;}
            let mut mn=DVec3::splat(f64::INFINITY);let mut mx=DVec3::splat(f64::NEG_INFINITY);
            for &vi in &self.shapes[i].sub_shapes{
                if let(Some(bmin),Some(bmax))=(self.shapes[vi].box_min,self.shapes[vi].box_max){mn=mn.min(bmin);mx=mx.max(bmax);}
            }
            if mn.x.is_finite(){self.shapes[i].box_min=Some(mn);self.shapes[i].box_max=Some(mx);self.shapes[i].box_gap=tol;}
        }
    }
    fn prepare_faces(&mut self,tol:f64){
        for fi in 0..self.nb_source_shapes{
            if self.shapes[fi].shape_type!=ShapeType::Face{continue;}
            let mut ns:HashSet<usize>=HashSet::new();
            for &wi in &self.shapes[fi].sub_shapes.clone(){
                if wi>=self.nb_shapes(){continue;}
                for &ei in &self.shapes[wi].sub_shapes.clone(){
                    if ei>=self.nb_shapes(){continue;}
                    if self.shapes[ei].shape_type==ShapeType::Edge{ns.insert(ei);for &vi in &self.shapes[ei].sub_shapes{if vi<self.nb_shapes(){ns.insert(vi);}}}
                }
            }
            self.shapes[fi].sub_shapes=ns.into_iter().collect();
            let mut mn=DVec3::splat(f64::INFINITY);let mut mx=DVec3::splat(f64::NEG_INFINITY);
            for &ss in &self.shapes[fi].sub_shapes{if let(Some(bmin),Some(bmax))=(self.shapes[ss].box_min,self.shapes[ss].box_max){mn=mn.min(bmin);mx=mx.max(bmax);}}
            if mn.x.is_finite(){self.shapes[fi].box_min=Some(mn);self.shapes[fi].box_max=Some(mx);self.shapes[fi].box_gap+=tol;}
        }
    }
    fn prepare_solids(&mut self){
        if self.arguments.len()!=1{return;}
        for si in 0..self.nb_source_shapes{
            if self.shapes[si].shape_type!=ShapeType::Solid{continue;}
            let mut ns:HashSet<usize>=HashSet::new();
            for &shi in &self.shapes[si].sub_shapes.clone(){if shi>=self.nb_shapes(){continue;}if self.shapes[shi].shape_type!=ShapeType::Shell{continue;}
                for &fi in &self.shapes[shi].sub_shapes{if fi>=self.nb_shapes(){continue;}ns.insert(fi);for &ei in &self.shapes[fi].sub_shapes{ns.insert(ei);}}}
            self.shapes[si].sub_shapes=ns.into_iter().collect();
        }
    }
    fn build_vertex_edge_map(&mut self){
        for ei in 0..self.nb_source_shapes{if self.shapes[ei].shape_type!=ShapeType::Edge{continue;}
            for &vi in &self.shapes[ei].sub_shapes{if vi>=self.nb_shapes(){continue;}let e=self.map_ve.entry(vi).or_default();if!e.contains(&ei){e.push(ei);}}}
    }
    fn build_bnd_box(&mut self,i:usize)->Option<(DVec3,DVec3,f64)>{
        if let(Some(mn),Some(mx))=(self.shapes[i].box_min,self.shapes[i].box_max){return Some((mn,mx,self.shapes[i].box_gap));}
        match self.shapes[i].shape_type{
            ShapeType::Vertex=>{
                let p=self.vertex_point(self.shapes[i].shape);
                let t=self.vertex_tolerance(self.shapes[i].shape);
                if let Some(pt)=p{let tol=t.max(1e-10);let b=(pt-DVec3::splat(tol),pt+DVec3::splat(tol),tol);self.shapes[i].box_min=Some(b.0);self.shapes[i].box_max=Some(b.1);self.shapes[i].box_gap=b.2;Some(b)}else{None}
            }
            _=>{
                let mut mn=DVec3::splat(f64::INFINITY);let mut mx=DVec3::splat(f64::NEG_INFINITY);let mut gap=0.0f64;
                for &c in &self.shapes[i].sub_shapes.clone(){if c<self.nb_shapes(){if let Some(b)=self.build_bnd_box(c){mn=mn.min(b.0);mx=mx.max(b.1);gap=gap.max(b.2);}}}
                if mn.x.is_finite(){self.shapes[i].box_min=Some(mn);self.shapes[i].box_max=Some(mx);self.shapes[i].box_gap=gap;Some((mn,mx,gap))}else{None}
            }
        }
    }
    fn vertex_tolerance(&self,s:Shape)->f64{s.vertex(&self.store).map_or(0.0,|vd|vd.tolerance)}
    fn vertex_point(&self,s:Shape)->Option<DVec3>{s.vertex(&self.store).map(|vd|vd.point)}

    // per-type count helpers
    pub fn vertex_count(&self)->usize{self.shapes.iter().filter(|s|s.shape_type==ShapeType::Vertex).count()}
    pub fn edge_count(&self)->usize{self.shapes.iter().filter(|s|s.shape_type==ShapeType::Edge).count()}
    pub fn face_count(&self)->usize{self.shapes.iter().filter(|s|s.shape_type==ShapeType::Face).count()}
    pub fn a_vertex_count(&self)->usize{self.shapes[..self.nb_source_shapes].iter().filter(|s|s.shape_type==ShapeType::Vertex).count()}
    pub fn a_edge_count(&self)->usize{self.shapes[..self.nb_source_shapes].iter().filter(|s|s.shape_type==ShapeType::Edge).count()}
    pub fn a_face_count(&self)->usize{self.shapes[..self.nb_source_shapes].iter().filter(|s|s.shape_type==ShapeType::Face).count()}
}

impl Default for DS{fn default()->Self{Self::new()}}
