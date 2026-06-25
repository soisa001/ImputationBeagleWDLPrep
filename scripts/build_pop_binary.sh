#!/usr/bin/env bash
# Build the static-musl pop-glimpse2 binary (no sudo, no docker). The Cromwell
# task VM can't reach crates.io inside the VPC-SC perimeter, so we build once on
# the notebook VM and pass the prebuilt binary via POP_BINARY_LOCAL.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # repo root
RS="${HERE}/pop_glimpse2_rust/pop-glimpse2.rs"
CARGO="${HERE}/pop_glimpse2_rust/Cargo.toml"
BUILD="${HOME}/pop-build"

[ -s "$RS" ] && [ -s "$CARGO" ] || { echo "ERROR: pop source missing under ${HERE}/pop_glimpse2_rust/"; exit 1; }

# rustc toolchain (user-local)
if ! command -v cargo >/dev/null 2>&1; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
fi
# shellcheck disable=SC1091
source "$HOME/.cargo/env"
rustup target add x86_64-unknown-linux-musl

rm -rf "$BUILD"; mkdir -p "$BUILD/src/bin" "$BUILD/.cargo"
cp "$RS"    "$BUILD/src/bin/pop-glimpse2.rs"
cp "$CARGO" "$BUILD/Cargo.toml"
printf '\n[workspace]\n' >> "$BUILD/Cargo.toml"          # don't climb to a stray ~/Cargo.toml
printf '[target.x86_64-unknown-linux-musl]\nlinker = "rust-lld"\n' > "$BUILD/.cargo/config.toml"

( cd "$BUILD" && rm -f Cargo.lock && cargo build --release --target x86_64-unknown-linux-musl )

BIN="$BUILD/target/x86_64-unknown-linux-musl/release/pop-glimpse2"
[ -x "$BIN" ] || { echo "ERROR: build produced no binary"; exit 1; }
echo "built: $BIN"
ldd "$BIN" 2>&1 | head -1 || true        # expect: not a dynamic executable
echo
echo "export POP_BINARY_LOCAL=\"$BIN\""
