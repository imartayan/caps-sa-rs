fn main() {
    println!("cargo:rerun-if-changed=cxx/CMakeLists.txt");
    println!("cargo:rerun-if-changed=cxx/wrapper.cpp");

    let mut cfg = cmake::Config::new("cxx");
    cfg.always_configure(false);

    if let Ok(cc_var) = std::env::var("CC") {
        cfg.define("CMAKE_C_COMPILER", cc_var);
    }
    if let Ok(cxx_var) = std::env::var("CXX") {
        cfg.define("CMAKE_CXX_COMPILER", cxx_var);
    }
    if cfg!(feature = "quiet") {
        cfg.define("CAPSSA_QUIET", "TRUE");
    }

    let dst = cfg.build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=caps_sa_wrapper");
    println!("cargo:rustc-link-lib=static=core");

    #[cfg(target_os = "linux")]
    println!("cargo:rustc-link-lib=dylib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=dylib=c++");
}
