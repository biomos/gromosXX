use std::env;
use std::process::Command;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    // Generate C header file using cbindgen (with config file)
    // Note: Binary I/O modules are Rust-only and excluded from C bindings
    let config = cbindgen::Config::from_file("cbindgen.toml").unwrap_or_default();
    match cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate()
    {
        Ok(bindings) => {
            bindings.write_to_file("gromos_rs.h");
            println!("cargo:warning=C bindings generated successfully");
        }
        Err(e) => {
            println!("cargo:warning=Unable to generate C bindings: {:?}", e);
            println!("cargo:warning=Binary I/O modules are Rust-only and excluded from C API");
        }
    }

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
