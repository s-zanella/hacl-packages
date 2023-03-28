#[cfg(not(windows))]
extern crate bindgen;

use std::{env, path::Path, process::Command};

#[cfg(not(windows))]
fn create_bindings(include_path: &Path, home_dir: &Path) {
    // Include paths
    let hacl_includes = vec![
        format!("-I{}", include_path.display()),
        format!("-I{}", include_path.join("hacl").display()),
        format!("-I{}", include_path.join("krml").display()),
        format!("-I{}", include_path.join("vale").display()),
    ];

    let bindings = bindgen::Builder::default()
        // Header to wrap HACL/Evercrypt headers
        .header("wrapper.h")
        // Set include paths for HACL/Evercrypt headers
        .clang_args(hacl_includes.iter())
        // Allow function we want to have in
        .allowlist_function("EverCrypt_AutoConfig2_.*")
        .allowlist_function("EverCrypt_AEAD_.*")
        .allowlist_function("EverCrypt_Curve25519_.*")
        .allowlist_function("EverCrypt_Ed25519_.*")
        .allowlist_function("EverCrypt_Hash_.*")
        .allowlist_function("EverCrypt_HKDF_.*")
        .allowlist_function("EverCrypt_HMAC_.*")
        .allowlist_function("Hacl_P256_.*")
        .allowlist_function("Hacl_RSAPSS_.*")
        .allowlist_function("Hacl_SHA3_.*")
        .allowlist_function("Hacl_Chacha20Poly1305_.*")
        .allowlist_function("Hacl_Hash_.*")
        .allowlist_function("Hacl_Streaming_.*")
        .allowlist_function("Hacl_Hash_Blake2.*")
        .allowlist_function("Hacl_Curve25519_.*")
        .allowlist_function("Hacl_HKDF_.*")
        .allowlist_function("Hacl_HMAC_.*")
        .allowlist_function("Hacl_HMAC_DRBG_.*")
        .allowlist_function("Hacl_Bignum64_.*")
        .allowlist_function("Hacl_Ed25519_.*")
        .allowlist_var("EverCrypt_Error_.*")
        .allowlist_var("Spec_.*")
        .allowlist_type("Spec_.*")
        .allowlist_type("Hacl_Streaming_SHA2_state.*")
        .allowlist_type("Hacl_Streaming_SHA3_state.*")
        .allowlist_type("Hacl_HMAC_DRBG_.*")
        // Block everything we don't need or define ourselves.
        .blocklist_type("EverCrypt_AEAD_state_s.*")
        // These functions currently use FFI-unsafe u128
        .blocklist_type("FStar_UInt128_uint128")
        .blocklist_function("Hacl_Hash_SHA2_update_last_384")
        .blocklist_function("Hacl_Hash_SHA2_update_last_512")
        .blocklist_function("Hacl_Hash_Blake2b_update_multi")
        .blocklist_function("Hacl_Hash_Blake2b_update_last")
        .blocklist_function("Hacl_Hash_Blake2b_Simd256_update_multi")
        .blocklist_function("Hacl_Hash_Blake2b_Simd256_update_last")
        // Disable tests to avoid warnings and keep it portable
        .layout_tests(false)
        // Generate bindings
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate bindings");

    // let bindings_path = out_path.join("bindings.rs");
    let home_bindings = home_dir.join("src/bindings/bindings.rs");
    bindings
        .write_to_file(home_bindings)
        .expect("Couldn't write bindings!");
}

#[cfg(windows)]
fn create_bindings(_: &Path, _: &Path) {}

fn build_hacl_c(path: &Path, cross_target: Option<String>) {
    println!(" >>> Building HACL C in {}", path.display());
    // cmake
    let mut cmake_cmd = Command::new("cmake");

    // Map cross compile targets to cmake toolchain files
    let toolchain_file = cross_target
        .map(|s| match s.as_str() {
            "x86_64-apple-darwin" => "-DCMAKE_TOOLCHAIN_FILE=config/x64-darwin.cmake",
            "aarch64-apple-darwin" => "-DCMAKE_TOOLCHAIN_FILE=config/aarch64-darwin.cmake",
            _ => "",
        })
        .unwrap_or_default();

    // We always build the release version here.
    // TODO: For debugging don't use this.
    let cmake_status = cmake_cmd
        .current_dir(path)
        .args(&[
            "-B",
            "build",
            "-G",
            "Ninja",
            "-D",
            "CMAKE_BUILD_TYPE=Release",
            toolchain_file,
        ])
        .status()
        .expect("Failed to run cmake.");
    if !cmake_status.success() {
        panic!("Failed to run cmake.")
    }
    // build
    let mut ninja_cmd = Command::new("ninja");
    let ninja_status = ninja_cmd
        .current_dir(path)
        .args(&["-f", "build.ninja", "-C", "build"])
        .status()
        .expect("Failed to run ninja.");
    if !ninja_status.success() {
        panic!("Failed to run ninja.")
    }

    // install
    let install_path = path.join("build").join("installed");
    println!(" >>> Installing HACL C into {}", install_path.display());
    let mut cmake_cmd = Command::new("cmake");
    let cmake_status = cmake_cmd
        .current_dir(path)
        .args(&[
            "--install",
            "build",
            "--prefix",
            install_path.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to install C library.");
    if !cmake_status.success() {
        panic!("Failed to install C library.")
    }
}

#[cfg(not(windows))]
fn copy_hacl_to_out(out_dir: &Path) {
    let mkdir_status = Command::new("mkdir")
        .arg("-p")
        .arg(out_dir.join("build"))
        .status()
        .expect("Failed to create build dir in out_dir.");
    if !mkdir_status.success() {
        panic!("Failed to create build dir in out_dir.")
    }

    fn copy_dir(path: &Path, out_dir: &Path) {
        let cp_status = Command::new("cp")
            .arg("-r")
            .arg(path)
            .arg(out_dir)
            .status()
            .expect("Failed to copy hacl to out_dir.");
        if !cp_status.success() {
            panic!("Failed to copy hacl to out_dir.")
        }
    }
    let local_c_path = Path::new(".c");
    copy_dir(&local_c_path.join("config"), &out_dir.join("config"));
    copy_dir(&local_c_path.join("src"), &out_dir.join("src"));
    copy_dir(&local_c_path.join("vale"), &out_dir.join("vale"));
    copy_dir(&local_c_path.join("karamel"), &out_dir.join("karamel"));
    copy_dir(&local_c_path.join("include"), &out_dir.join("include"));
    copy_dir(
        &local_c_path.join("config").join("default_config.cmake"),
        &out_dir.join("build").join("config.cmake"),
    );
    copy_dir(&local_c_path.join("CMakeLists.txt"), out_dir);
}

// TODO: make this work on windows.
#[cfg(windows)]
fn copy_hacl_to_out(out_dir: &Path) {
    panic!("TODO: Windows build is not supported right now.");
    // let cp_status = Command::new("cmd")
    //     .args(&[
    //         "/C",
    //         "robocopy",
    //         ".c",
    //         &format!("{}\\.c", out_dir.to_str().unwrap()),
    //         "/e",
    //         "/s",
    //     ])
    //     .status()
    //     .expect(&format!("Failed to copy hacl to {:?}", out_dir));

    // println!("Return code {}", cp_status.code().unwrap());

    // println!("Copied hacl-star to {:?}", out_dir);
}

fn main() {
    // Get ENV variables
    let home_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let home_dir = Path::new(&home_dir);
    let mach_build = env::var("MACH_BUILD").ok().is_some();
    let target = env::var("TARGET").unwrap();
    let host = env::var("HOST").unwrap();
    let out_dir = env::var("OUT_DIR").unwrap();
    let out_dir = Path::new(&out_dir);
    println!("mach_build: {}", mach_build);

    let cross_target = if target != host { Some(target) } else { None };

    // Get the C library and build it first.
    // This is the default behaviour. It can be disabled when working on this
    // to pick up the local version. This is what the global mach script does.
    let hacl_path = if !mach_build {
        // Copy all of the code into out to prepare build
        let c_out_dir = out_dir.join("c");
        println!(" >>> Copying HACL C file");
        println!("     from {}", home_dir.join(".c").display());
        println!("     to {}", c_out_dir.display());
        copy_hacl_to_out(&c_out_dir);
        build_hacl_c(&c_out_dir, cross_target);

        c_out_dir.join("build").join("installed")
    } else {
        // Use the higher level install directory.
        home_dir
            .join("..")
            .join("..")
            .join("build")
            .join("installed")
    };
    let hacl_lib_path = hacl_path.join("lib");
    let hacl_include_path = hacl_path.join("include");

    // Set library name to look up
    let library_name = "hacl_static";

    // Set re-run trigger
    println!("cargo:rerun-if-changed=wrapper.h");
    // We should re-run if the library changed. But this triggers the build
    // to re-run every time right now.
    // println!(
    //     "cargo:rerun-if-changed={}",
    //     hacl_lib_path.join(library_name).display()
    // );

    // Generate new bindings. This is a no-op on Windows.
    create_bindings(&hacl_include_path, home_dir);

    // Link hacl library.
    let mode = "static";
    println!("cargo:rustc-link-lib={}={}", mode, library_name);
    println!("cargo:rustc-link-search=native={}", hacl_lib_path.display());
    println!("cargo:lib={}", hacl_lib_path.display());
}