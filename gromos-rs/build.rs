use std::env;

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

    println!("cargo:rerun-if-changed=src/");
}
