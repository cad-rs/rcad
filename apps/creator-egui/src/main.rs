#[cfg(not(target_arch = "wasm32"))]
fn print_usage(program: &str) {
    eprintln!("Usage: {program} [--step <file.step> | <file.step>]");
}

#[cfg(not(target_arch = "wasm32"))]
fn parse_step_arg(args: &[String]) -> Result<Option<String>, String> {
    let mut positional_path: Option<&str> = None;
    let mut iter = args.iter().skip(1);

    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--help" | "-h" => {
                return Err(String::new());
            }
            "--step" | "-s" => {
                let Some(path) = iter.next() else {
                    return Err("missing file path after --step".to_string());
                };
                if positional_path.is_some() {
                    return Err("STEP file specified more than once".to_string());
                }
                positional_path = Some(path);
            }
            value if value.starts_with('-') => {
                return Err(format!("unknown option: {value}"));
            }
            value => {
                if positional_path.is_some() {
                    return Err("STEP file specified more than once".to_string());
                }
                positional_path = Some(value);
            }
        }
    }

    positional_path
        .map(|path| {
            std::fs::read_to_string(path)
                .map_err(|err| format!("failed to read STEP file '{path}': {err}"))
        })
        .transpose()
}

fn main() {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let args: Vec<String> = std::env::args().collect();
        let program = args
            .first()
            .map(String::as_str)
            .unwrap_or("creator-egui");

        let step_content = match parse_step_arg(&args) {
            Ok(content) => content,
            Err(message) if message.is_empty() => {
                print_usage(program);
                std::process::exit(0);
            }
            Err(message) => {
                eprintln!("{message}");
                print_usage(program);
                std::process::exit(2);
            }
        };

        creator_egui::run_native(step_content);
    }
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::parse_step_arg;

    #[test]
    fn parse_step_arg_accepts_positional_path() {
        let args = vec!["creator-egui".to_string(), "assets/box.step".to_string()];
        let content = parse_step_arg(&args).expect("positional path should parse");
        assert!(content.is_some());
    }

    #[test]
    fn parse_step_arg_accepts_flagged_path() {
        let args = vec![
            "creator-egui".to_string(),
            "--step".to_string(),
            "assets/box.step".to_string(),
        ];
        let content = parse_step_arg(&args).expect("--step path should parse");
        assert!(content.is_some());
    }

    #[test]
    fn parse_step_arg_rejects_duplicate_paths() {
        let args = vec![
            "creator-egui".to_string(),
            "--step".to_string(),
            "assets/box.step".to_string(),
            "assets/box.step".to_string(),
        ];
        let err = parse_step_arg(&args).expect_err("duplicate paths should fail");
        assert!(err.contains("more than once"));
    }
}
