# caps-sa-rs

Rust bindings for [CaPS-SA](https://github.com/jamshed/CaPS-SA), a parallel and cache-friendly Suffix Array and LCP array construction algorithm.

[Documentation](https://imartayan.github.io/caps-sa-rs/)

## Usage

```toml
[dependencies]
caps-sa-rs = { git = "https://github.com/imartayan/caps-sa-rs.git" }
```

When using external memory on macOS, you might have to increase the limit of opened files:
```sh
ulimit -n 2048
```
