use std::env;
use std::process::Command;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    // Generate C header file using cbindgen
    cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_language(cbindgen::Language::C)
        .with_include_guard("GROMOS_RS_H")
        .with_documentation(true)
        .with_pragma_once(true)
        .generate()
        .expect("Unable to generate C bindings")
        .write_to_file("gromos_rs.h");

    // Auto-detect FFTW3 library and enable feature if available
    if is_fftw3_available() {
        println!("cargo:rustc-cfg=feature=\"use-fftw\"");
        println!("cargo:warning=FFTW3 detected - using high-performance FFT");
    } else {
        println!("cargo:warning=FFTW3 not found - using rustfft fallback");
    }

    println!("cargo:rerun-if-changed=src/");
}

/// Check if FFTW3 library is available on the system
fn is_fftw3_available() -> bool {
    // Try pkg-config first (Linux, macOS with pkg-config)
    if Command::new("pkg-config")
        .args(&["--exists", "fftw3"])
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
    {
        return true;
    }

    // Try to find library manually (for systems without pkg-config)
    // Check common library paths
    let lib_paths = [
        "/usr/lib",
        "/usr/local/lib",
        "/opt/homebrew/lib",  // macOS ARM
        "/usr/lib/x86_64-linux-gnu",  // Debian/Ubuntu
    ];

    for path in &lib_paths {
        let lib_file_so = format!("{}/libfftw3.so", path);
        let lib_file_dylib = format!("{}/libfftw3.dylib", path);
        let lib_file_a = format!("{}/libfftw3.a", path);

        if std::path::Path::new(&lib_file_so).exists()
            || std::path::Path::new(&lib_file_dylib).exists()
            || std::path::Path::new(&lib_file_a).exists()
        {
            return true;
        }
    }

    false
}
