//! TKDESTEP-aligned tests for `clean_text_for_send` (OCCT `StepData_StepWriter::CleanTextForSend`).
//!
//! OCCT source: src/DataExchange/TKDESTEP/GTests/StepData_StepWriter_Test.cxx

use rcad_step::clean_text_for_send;

// ── Basic escaping ──

#[test]
fn cleans_single_quotes() {
    assert_eq!(
        clean_text_for_send("text with 'single quotes'"),
        "text with ''single quotes''"
    );
}

#[test]
fn cleans_backslashes() {
    assert_eq!(
        clean_text_for_send("path\\with\\backslashes"),
        "path\\\\with\\\\backslashes"
    );
}

#[test]
fn cleans_newlines() {
    assert_eq!(clean_text_for_send("line1\nline2"), "line1\\N\\line2");
}

#[test]
fn cleans_tabs() {
    assert_eq!(
        clean_text_for_send("text\twith\ttabs"),
        "text\\T\\with\\T\\tabs"
    );
}

// ── Control directive preservation ──

#[test]
fn preserves_x_directive() {
    assert_eq!(
        clean_text_for_send("text with \\XA7\\ section sign"),
        "text with \\XA7\\ section sign"
    );
}

#[test]
fn preserves_x2_directive() {
    assert_eq!(
        clean_text_for_send("\\X2\\03C0\\X0\\ is pi"),
        "\\X2\\03C0\\X0\\ is pi"
    );
}

#[test]
fn preserves_x4_directive() {
    assert_eq!(
        clean_text_for_send("emoji \\X4\\001F600\\X0\\ face"),
        "emoji \\X4\\001F600\\X0\\ face"
    );
}

#[test]
fn preserves_s_directive() {
    assert_eq!(
        clean_text_for_send("text with \\S\\ directive"),
        "text with \\S\\ directive"
    );
}

#[test]
fn preserves_p_directive() {
    assert_eq!(
        clean_text_for_send("\\PA\\ code page setting"),
        "\\PA\\ code page setting"
    );
}

// ── Existing directive preservation ──

#[test]
fn preserves_existing_n_directive() {
    assert_eq!(clean_text_for_send("line1\\N\\line2"), "line1\\N\\line2");
}

#[test]
fn preserves_existing_t_directive() {
    assert_eq!(
        clean_text_for_send("text\\T\\with\\T\\tab"),
        "text\\T\\with\\T\\tab"
    );
}

// ── Mixed content ──

#[test]
fn mixed_quotes_and_directives() {
    assert_eq!(
        clean_text_for_send("see \\XA7\\ section and 'quotes'"),
        "see \\XA7\\ section and ''quotes''"
    );
}

#[test]
fn mixed_backslashes_and_directives() {
    assert_eq!(
        clean_text_for_send("\\XA7\\ and path\\file"),
        "\\XA7\\ and path\\\\file"
    );
}

#[test]
fn mixed_directive_quotes_tab() {
    assert_eq!(
        clean_text_for_send("prefix \\X2\\03B103B2\\X0\\ 'text' with\ttab"),
        "prefix \\X2\\03B103B2\\X0\\ ''text'' with\\T\\tab"
    );
}

// ── Edge cases ──

#[test]
fn empty_string() {
    assert_eq!(clean_text_for_send(""), "");
}

#[test]
fn only_quotes() {
    assert_eq!(clean_text_for_send("''"), "''''");
}

#[test]
fn only_control_directive() {
    assert_eq!(clean_text_for_send("\\XA7\\"), "\\XA7\\");
}

#[test]
fn consecutive_directives() {
    assert_eq!(clean_text_for_send("\\XA7\\\\XB6\\"), "\\XA7\\\\XB6\\");
}

// ── Malformed but safe input ──

#[test]
fn incomplete_directive_with_quotes() {
    assert_eq!(
        clean_text_for_send("incomplete \\X and 'quotes'"),
        "incomplete \\\\X and ''quotes''"
    );
}

#[test]
fn partial_directive() {
    assert_eq!(
        clean_text_for_send("partial \\XA and more"),
        "partial \\\\XA and more"
    );
}

// ── Hex sequence detection ──

#[test]
fn x2_hex_sequence() {
    assert_eq!(
        clean_text_for_send("\\X2\\03B103B203B3\\X0\\"),
        "\\X2\\03B103B203B3\\X0\\"
    );
}

#[test]
fn x4_hex_sequence() {
    assert_eq!(
        clean_text_for_send("\\X4\\001F600001F638\\X0\\"),
        "\\X4\\001F600001F638\\X0\\"
    );
}

#[test]
fn hex_sequence_with_surrounding_text() {
    assert_eq!(
        clean_text_for_send("start \\X2\\03C0\\X0\\ end"),
        "start \\X2\\03C0\\X0\\ end"
    );
}
