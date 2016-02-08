extern crate syntex;
extern crate serde_codegen;

use std::env;
use std::path::Path;

pub fn main() {
    generate_serde_traits("src/data_structures/fmindex.rs.in", "fmindex.rs");
    generate_serde_traits("src/data_structures/bwt.rs.in", "bwt.rs");
}

fn generate_serde_traits(source_file: &'static str, dest_file: &'static str) {
    let out_dir = env::var_os("OUT_DIR").unwrap();

    let src = Path::new(source_file);
    let dst = Path::new(&out_dir).join(dest_file);

    let mut registry = syntex::Registry::new();
    serde_codegen::register(&mut registry);
    registry.expand("", &src, &dst).unwrap();
}
