set -ve
prefix=bam2msa
base=~/webassembly_crystal/
$base/crystal/bin/crystal build $prefix.cr --cross-compile --target wasm32-unknown-wasi --release
/usr/lib/llvm-12/bin/wasm-ld  $prefix.wasm  $base/wasi-sysroot/lib/wasm32-wasi/crt1.o  -o ${prefix}_final.wasm -L $base/wasi-sysroot/lib/wasm32-wasi -lc -lpcre -lclang_rt.builtins-wasm32 
