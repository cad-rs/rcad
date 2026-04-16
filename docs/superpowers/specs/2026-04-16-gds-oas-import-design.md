# GDS/OAS文件导入导出功能设计

**日期:** 2026-04-16
**状态:** 设计中
**依赖:** laykit crate

## 1. 概述

### 1.1 目标

为RCAD CAD引擎添加GDSII (GDS) 和 OASIS (OAS) 文件格式的导入导出功能，支持：
- 将2D版图数据转换为3D实体模型
- 通过配置文件指定各层厚度
- 保持层级结构（cell引用关系）
- 支持导出回GDS/OAS格式

### 1.2 使用场景

- VLSI/IC设计的3D可视化
- MEMS结构建模
- 光子集成电路设计
- 生成用于仿真/制造的3D实体模型

### 1.3 laykit能力

laykit crate提供完整的GDS/OAS支持：
- **GDSII**: 完整读写支持，支持全部7种元素类型
- **OASIS**: 完整读写支持
- **布尔运算**: union, intersection, difference, xor, slice, offset
- **拓扑处理**: cell层次结构、展开、依赖分析
- **流式解析**: 支持大文件（>1GB）

## 2. 需求规格

| 需求项 | 决定 |
|--------|------|
| 用途 | 生成3D实体用于仿真/制造 |
| 厚度处理 | 通过配置文件/参数指定每层厚度 |
| 格式支持 | GDS + OAS（laykit完整支持） |
| 层级处理 | 保持层级结构 |
| 模块位置 | 先建rcad-gds和rcad-oas，后统一到rcad-io |
| 实现范围 | 全功能：导入+拉伸+导出 |

## 3. 架构设计

### 3.1 整体架构

采用分步实现策略：

```
阶段1: rcad-gds crate  ─┐
                        ├─► 阶段3: rcad-io crate (统一接口)
阶段2: rcad-oas crate ─┘
```

### 3.2 模块依赖

```
rcad-gds / rcad-oas
    │
    ├─► laykit (GDS/OAS解析)
    │
    └─► rcad-kernel (BRep数据结构)
```

## 4. 数据结构设计

### 4.1 层配置

```rust
/// 层配置：指定每层的拉伸参数
pub struct LayerConfig {
    pub layers: HashMap<i32, LayerSettings>,
}

pub struct LayerSettings {
    pub thickness: f64,        // 厚度（用户单位）
    pub z_offset: f64,         // Z方向偏移
    pub color: Option<[f32; 4]>, // RGBA颜色
    pub name: Option<String>,   // 可选名称
}

impl LayerConfig {
    /// 从JSON文件加载
    pub fn from_json(path: &Path) -> Result<Self, ConfigError>;

    /// 创建默认配置
    pub fn default() -> Self;

    /// 添加层设置
    pub fn with_layer(mut self, layer: i32, settings: LayerSettings) -> Self;
}
```

### 4.2 GDS中间表示

```rust
/// GDS库（对应GDSII的Library）
pub struct GdsLibrary {
    pub name: String,
    pub units: GdsUnits,
    pub structures: HashMap<String, GdsStructure>,
}

pub struct GdsUnits {
    pub user_unit: f64,    // 用户单位
    pub meter_unit: f64,   // 米单位
}

pub struct GdsStructure {
    pub name: String,
    pub boundaries: Vec<GdsBoundary>,   // 多边形
    pub paths: Vec<GdsPath>,            // 路径
    pub texts: Vec<GdsText>,            // 文本标注
    pub references: Vec<GdsReference>,  // cell引用
}

pub struct GdsBoundary {
    pub layer: i32,
    pub datatype: i32,
    pub points: Vec<DVec2>,  // 闭合多边形顶点
}

pub struct GdsPath {
    pub layer: i32,
    pub datatype: i32,
    pub width: f64,
    pub points: Vec<DVec2>,  // 路径中心线
    pub end_cap: EndCapType,
}

pub struct GdsText {
    pub layer: i32,
    pub text_type: i32,
    pub position: DVec2,
    pub content: String,
}

pub struct GdsReference {
    pub cell_name: String,
    pub transform: Transform2D,
    pub array: Option<ArrayParams>,
}

pub struct ArrayParams {
    pub columns: u32,
    pub rows: u32,
    pub column_offset: DVec2,
    pub row_offset: DVec2,
}
```

### 4.3 转换结果

```rust
/// 导入结果：保持层级结构
pub struct LayoutResult {
    pub top_cell: String,
    pub cells: HashMap<String, Cell3D>,
}

pub struct Cell3D {
    pub name: String,
    pub shapes: Vec<Shape3D>,        // 本cell的形状
    pub references: Vec<CellRef>,    // 对其他cell的引用
    pub bbox: BoundingBox,           // 边界框
}

pub struct Shape3D {
    pub layer: i32,
    pub brep: BRep,                  // 3D实体
    pub transform: Option<Transform>, // 变换
}

pub struct CellRef {
    pub cell_name: String,
    pub transform: Transform,        // 位置和旋转
    pub array: Option<ArrayParams>,  // 是否为阵列引用
}
```

## 5. API设计

### 5.1 GDS读取API

```rust
pub struct GdsReader;

impl GdsReader {
    /// 从文件读取GDS
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<GdsLibrary, GdsError>;

    /// 从字节数据解析
    pub fn parse_bytes(data: &[u8]) -> Result<GdsLibrary, GdsError>;
}

impl GdsLibrary {
    /// 获取顶层cell（未被其他cell引用的cell）
    pub fn top_cells(&self) -> Vec<&str>;

    /// 转换为BRep（扁平化，带层配置）
    pub fn to_brep(&self, cell_name: &str, config: &LayerConfig) -> Result<BRep, GdsError>;

    /// 转换为LayoutResult（保持层级结构）
    pub fn to_layout(&self, cell_name: &str, config: &LayerConfig) -> Result<LayoutResult, GdsError>;
}
```

### 5.2 GDS写入API

```rust
pub struct GdsWriter;

impl GdsWriter {
    /// 写入文件
    pub fn write_file<P: AsRef<Path>>(library: &GdsLibrary, path: P) -> Result<(), GdsError>;

    /// 转换为字节数据
    pub fn to_bytes(library: &GdsLibrary) -> Result<Vec<u8>, GdsError>;
}

impl GdsLibrary {
    /// 从BRep创建GDS库（反向转换）
    pub fn from_brep(brep: &BRep, layer_mapping: &LayerMapping) -> Self;
}
```

### 5.3 OASIS API

OASIS API结构与GDS类似：

```rust
pub struct OasReader;
impl OasReader {
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<OasLibrary, OasError>;
    pub fn parse_bytes(data: &[u8]) -> Result<OasLibrary, OasError>;
}

pub struct OasWriter;
impl OasWriter {
    pub fn write_file<P: AsRef<Path>>(library: &OasLibrary, path: P) -> Result<(), OasError>;
    pub fn to_bytes(library: &OasLibrary) -> Result<Vec<u8>, OasError>;
}
```

### 5.4 统一IO接口（rcad-io）

```rust
/// 统一的布局格式
pub enum LayoutFormat {
    Gds,
    Oasis,
}

pub struct LayoutIo;

impl LayoutIo {
    /// 自动检测格式并读取
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn LayoutLibrary>, IoError>;

    /// 写入指定格式
    pub fn write_file(library: &dyn LayoutLibrary, format: LayoutFormat, path: P) -> Result<(), IoError>;
}

/// 统一的布局库trait
pub trait LayoutLibrary {
    fn name(&self) -> &str;
    fn top_cells(&self) -> Vec<&str>;
    fn to_brep(&self, cell: &str, config: &LayerConfig) -> Result<BRep, IoError>;
}
```

## 6. 转换流程

### 6.1 GDS → BRep 转换流程

```
GDS文件
    │
    ▼ laykit解析
GdsLibrary (中间表示)
    │
    ├─ 选择顶层cell
    │
    ▼ 递归处理
┌──────────────────────────────┐
│  处理GdsStructure            │
│  ├─ Boundary → Wire → Face   │
│  ├─ Path → Wire → Face       │
│  ├─ Text → Annotation        │
│  └─ Reference → 递归处理     │
└──────────────────────────────┘
    │
    ▼ 应用LayerConfig
┌──────────────────────────────┐
│  2D Face → 3D Solid          │
│  ├─ 查询层厚度配置           │
│  ├─ 按厚度拉伸               │
│  └─ 应用Z偏移                │
└──────────────────────────────┘
    │
    ▼ 合并结果
BRep (3D实体)
```

### 6.2 BRep → GDS 转换流程

```
BRep (3D实体)
    │
    ▼ 分解为Faces
Vec<Face>
    │
    ▼ 投影到XY平面
Vec<Wire2D>
    │
    ├─ 按Z高度分组
    │
    ▼ 创建GdsStructure
┌──────────────────────────────┐
│  Wire → Boundary             │
│  ├─ 计算层号（从Z高度）      │
│  └─ 转换坐标                 │
└──────────────────────────────┘
    │
    ▼ laykit序列化
GDS文件
```

### 6.3 层级结构保持策略

```rust
impl GdsLibrary {
    pub fn to_layout(&self, cell_name: &str, config: &LayerConfig) -> Result<LayoutResult, GdsError> {
        let mut cells = HashMap::new();
        self.build_cell_recursive(cell_name, config, &mut cells)?;

        Ok(LayoutResult {
            top_cell: cell_name.to_string(),
            cells,
        })
    }

    fn build_cell_recursive(
        &self,
        name: &str,
        config: &LayerConfig,
        cells: &mut HashMap<String, Cell3D>
    ) -> Result<(), GdsError> {
        // 检查是否已处理（避免循环引用）
        if cells.contains_key(name) {
            return Ok(());
        }

        let structure = self.structures.get(name)
            .ok_or(GdsError::CellNotFound(name.to_string()))?;

        // 转换当前cell的几何
        let shapes = self.convert_shapes(structure, config)?;

        // 收集引用（不展开）
        let references: Vec<CellRef> = structure.references.iter()
            .map(|r| self.convert_reference(r))
            .collect();

        cells.insert(name.to_string(), Cell3D {
            name: name.to_string(),
            shapes,
            references,
            bbox: BoundingBox::default(),
        });

        // 递归处理被引用的cell
        for ref_ in &structure.references {
            self.build_cell_recursive(&ref_.cell_name, config, cells)?;
        }

        Ok(())
    }
}
```

## 7. 文件结构

### 7.1 rcad-gds crate

```
libs/rcad-gds/
├── Cargo.toml
└── src/
    ├── lib.rs              # 公开API入口
    ├── reader.rs           # GdsReader实现
    ├── writer.rs           # GdsWriter实现
    ├── convert.rs          # GdsLibrary <-> BRep 转换
    ├── types.rs            # GdsLibrary, GdsStructure等类型
    ├── layer_config.rs     # LayerConfig配置
    ├── error.rs            # GdsError错误类型
    └── geometry.rs         # 几何转换辅助函数
```

### 7.2 rcad-oas crate

```
libs/rcad-oas/
├── Cargo.toml
└── src/
    ├── lib.rs              # 公开API入口
    ├── reader.rs           # OasReader实现
    ├── writer.rs           # OasWriter实现
    ├── convert.rs          # OasLibrary <-> BRep 转换
    ├── types.rs            # OasLibrary等类型
    ├── layer_config.rs     # 复用或独立实现
    ├── error.rs            # OasError错误类型
    └── geometry.rs         # 几何转换辅助函数
```

### 7.3 rcad-io crate（第三阶段）

```
libs/rcad-io/
├── Cargo.toml
└── src/
    ├── lib.rs              # 统一API：LayoutIo
    ├── format.rs           # LayoutFormat枚举
    ├── traits.rs           # LayoutLibrary trait
    ├── detection.rs        # 格式自动检测
    └── layer_config.rs     # 统一LayerConfig（re-export）
```

## 8. 依赖配置

### 8.1 rcad-gds/Cargo.toml

```toml
[package]
name = "rcad-gds"
version = "0.1.0"
edition = "2024"

[dependencies]
laykit = "0.1"
rcad-kernel = { path = "../rcad-kernel" }
glam = { workspace = true }
thiserror = "2.0"
serde = { version = "1.0", features = ["derive"] }

[dev-dependencies]
tempfile = "3"
```

### 8.2 rcad-oas/Cargo.toml

```toml
[package]
name = "rcad-oas"
version = "0.1.0"
edition = "2024"

[dependencies]
laykit = "0.1"
rcad-kernel = { path = "../rcad-kernel" }
glam = { workspace = true }
thiserror = "2.0"
serde = { version = "1.0", features = ["derive"] }

[dev-dependencies]
tempfile = "3"
```

## 9. 错误处理

```rust
#[derive(Debug, thiserror::Error)]
pub enum GdsError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid GDS format: {0}")]
    InvalidFormat(String),

    #[error("Cell not found: {0}")]
    CellNotFound(String),

    #[error("Layer {0} not configured")]
    LayerNotConfigured(i32),

    #[error("Geometry conversion failed: {0}")]
    GeometryError(String),

    #[error("Empty structure: {0}")]
    EmptyStructure(String),

    #[error("laykit error: {0}")]
    Laykit(String),
}
```

## 10. 测试策略

### 10.1 单元测试

- `reader_test.rs`：GDS解析正确性
- `writer_test.rs`：GDS序列化正确性
- `convert_test.rs`：几何转换正确性
- `layer_config_test.rs`：配置解析正确性

### 10.2 集成测试

- `tests/integration/`
  - `read_write_roundtrip.rs`：读写往返测试
  - `gds_to_brep.rs`：GDS转BRep测试
  - `brep_to_gds.rs`：BRep转GDS测试

### 10.3 测试数据

- `tests/data/`
  - `simple_boundary.gds`：简单多边形
  - `hierarchical.gds`：层级结构
  - `with_paths.gds`：包含路径
  - `complex.gds`：复杂综合测试

## 11. 实现阶段

### 阶段1：rcad-gds（预计工作量：中）
1. 创建crate骨架
2. 实现GdsReader（基于laykit）
3. 实现LayerConfig
4. 实现GDS → BRep转换
5. 实现BRep → GDS转换
6. 实现GdsWriter
7. 编写测试

### 阶段2：rcad-oas（预计工作量：小）
1. 复制rcad-gds结构
2. 适配OASIS格式差异
3. 实现读写转换
4. 编写测试

### 阶段3：rcad-io（预计工作量：小）
1. 创建统一接口
2. 实现格式自动检测
3. 定义LayoutLibrary trait
4. 整合两个crate

## 12. 风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| laykit API变化 | 可能需要适配代码 | 锁定版本，关注更新 |
| 大文件内存占用 | 可能导致OOM | 使用laykit的流式解析 |
| 复杂几何转换失败 | 部分图形无法转换 | 提供错误详情，跳过并记录 |
| 循环引用 | 无限递归 | 已处理：检查是否已访问 |
