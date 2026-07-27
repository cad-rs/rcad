// rcad/libs/rcad-gds/tests/gds_roundtrip.rs

use glam::DVec2;
use rcad_gds::{
    GdsBoundary, GdsLibrary, GdsReader, GdsReference, GdsStructure, GdsWriter, LayerConfig,
    LayerSettings, Transform2D,
};

/// Test basic roundtrip: create library -> write -> read -> compare
#[test]
fn test_basic_roundtrip() {
    let mut library = GdsLibrary {
        name: "ROUNDTRIP_TEST".to_string(),
        ..Default::default()
    };

    let structure = GdsStructure {
        name: "TOP".to_string(),
        boundaries: vec![
            GdsBoundary {
                layer: 1,
                datatype: 0,
                points: vec![
                    DVec2::new(0.0, 0.0),
                    DVec2::new(1000.0, 0.0),
                    DVec2::new(1000.0, 1000.0),
                    DVec2::new(0.0, 1000.0),
                    DVec2::new(0.0, 0.0),
                ],
            },
            GdsBoundary {
                layer: 2,
                datatype: 0,
                points: vec![
                    DVec2::new(100.0, 100.0),
                    DVec2::new(900.0, 100.0),
                    DVec2::new(900.0, 900.0),
                    DVec2::new(100.0, 900.0),
                    DVec2::new(100.0, 100.0),
                ],
            },
        ],
        ..Default::default()
    };

    library.structures.insert("TOP".to_string(), structure);

    // Write to bytes
    let bytes = GdsWriter::to_bytes(&library).expect("Failed to write");

    // Read back
    let restored = GdsReader::parse_bytes(&bytes).expect("Failed to parse");

    // Verify
    assert_eq!(restored.name, library.name);
    assert!(restored.has_cell("TOP"));

    let top = restored.structures.get("TOP").unwrap();
    assert_eq!(top.boundaries.len(), 2);
}

/// Test conversion to BRep
#[test]
fn test_to_brep() {
    let mut library = GdsLibrary {
        name: "BREP_TEST".to_string(),
        ..Default::default()
    };

    let structure = GdsStructure {
        name: "TOP".to_string(),
        boundaries: vec![GdsBoundary {
            layer: 1,
            datatype: 0,
            points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(100.0, 0.0),
                DVec2::new(100.0, 100.0),
                DVec2::new(0.0, 100.0),
                DVec2::new(0.0, 0.0),
            ],
        }],
        ..Default::default()
    };

    library.structures.insert("TOP".to_string(), structure);

    let config = LayerConfig::new().with_layer(1, LayerSettings::new(10.0));

    let brep = library.to_brep("TOP", &config).expect("Failed to convert");

    assert!(!brep.solids.is_empty());
    assert!(!brep.vertices.is_empty());
}

/// Test hierarchical structure
#[test]
fn test_hierarchical() {
    let mut library = GdsLibrary {
        name: "HIER_TEST".to_string(),
        ..Default::default()
    };

    // Create leaf cell
    let leaf = GdsStructure {
        name: "LEAF".to_string(),
        boundaries: vec![GdsBoundary {
            layer: 1,
            datatype: 0,
            points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(10.0, 0.0),
                DVec2::new(10.0, 10.0),
                DVec2::new(0.0, 10.0),
                DVec2::new(0.0, 0.0),
            ],
        }],
        ..Default::default()
    };

    // Create top cell with reference to leaf
    let top = GdsStructure {
        name: "TOP".to_string(),
        references: vec![
            GdsReference {
                cell_name: "LEAF".to_string(),
                transform: Transform2D::from_translation(100.0, 0.0),
                array: None,
            },
            GdsReference {
                cell_name: "LEAF".to_string(),
                transform: Transform2D::from_translation(200.0, 0.0),
                array: None,
            },
        ],
        ..Default::default()
    };

    library.structures.insert("LEAF".to_string(), leaf);
    library.structures.insert("TOP".to_string(), top);

    // Verify top cells
    let top_cells = library.top_cells();
    assert_eq!(top_cells, vec!["TOP"]);

    // Roundtrip
    let bytes = GdsWriter::to_bytes(&library).expect("Failed to write");
    let restored = GdsReader::parse_bytes(&bytes).expect("Failed to parse");

    assert_eq!(restored.structures.len(), 2);

    let top_restored = restored.structures.get("TOP").unwrap();
    assert_eq!(top_restored.references.len(), 2);
}
