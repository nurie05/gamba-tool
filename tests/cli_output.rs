use std::process::Command;
use std::str;

#[test]
fn test_operon_summary_output() {
    // Run the binary (assumes the package provides a bin target, e.g. `gamba`)
    let output = Command::new(env!("CARGO_BIN_EXE_gamba"))
        .args([
            "-f", "tests/resources/Samples_test_chr1.gtf",
            "-o", ".tests/test_operon_summary_output",
        ])
        .output()
        .expect("Failed to execute gamba_tool binary");

    assert!(
        output.status.success(),
        "Program exited with error: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);

    // Helper closure to check if expected text appears (ignoring timestamps/log prefixes)
    let contains_line = |substr: &str| -> bool {
        stdout
            .lines()
            .any(|line| line.contains(substr))
    };

    // Assertions (ignoring timestamp/log level prefixes)
    assert!(contains_line("Processing chromosome Chr1 (11722 transcripts)"));
    assert!(contains_line("Total number of OPRNs found: 60"));
    assert!(contains_line("Total number of OpGs found: 125"));
    assert!(contains_line("2 genes: 55"));
    assert!(contains_line("3 genes: 5"));
    assert!(contains_line("4 genes: 0"));
    assert!(contains_line("5 genes: 0"));
    assert!(contains_line(">5 genes: 0"));
}
