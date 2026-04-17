# DNA Seed

Lossless data compressor that beats bzip2, xz, and zstd — in both size and speed.

Pure C, single file, no tricks. Just better algorithms.

## Benchmarks (100MB enwik8)

| Compressor | Size | Time | Savings |
|---|---|---|---|
| **DNA Seed v12** | **24,819,858** | **6.2s** | **75.2%** |
| xz -9 | 24,865,244 | 118s | 75.1% |
| zstd -19 | 26,954,633 | 83s | 73.0% |
| bzip2 -9 | 29,008,758 | 15s | 71.0% |
| gzip -9 | 36,445,248 | 7.6s | 63.6% |

Tested on Intel i7-9750H, 13GB RAM, Linux (WSL2).

## How it works

```
Input → Pre-BWT RLE → BWT (divsufsort) → MTF → ZRLE → Adaptive Arithmetic → Seed
```

- **BWT** (Burrows-Wheeler Transform) sorts all rotations of the input, grouping similar contexts together
- **MTF** (Move-to-Front) converts the BWT output into a stream dominated by zeros
- **ZRLE** (Zero Run-Length Encoding) compresses the zero runs using bijective base-2 coding
- **Adaptive Arithmetic Coding** with Krichevsky-Trofimov estimator and exponential decay — no frequency table needed, encoder and decoder learn the distribution as they go
- **Block-parallel** processing with pthreads — splits input into blocks, processes in parallel

The key innovation is the KT estimator with exponential decay (alpha=0.99961), implemented in fixed-point integer arithmetic. No floats in the hot path. This gives local adaptation comparable to bzip2's multi-table Huffman, with lower overhead.

## Build

Requires `libdivsufsort`:

```bash
# Debian/Ubuntu
sudo apt install libdivsufsort-dev

# macOS
brew install libdivsufsort
```

Then:

```bash
make
```

## Usage

```bash
# Compress
./dna_seed compress input.txt output.dna

# Decompress
./dna_seed decompress output.dna restored.txt

# Benchmark (compress + decompress + verify roundtrip)
./dna_seed benchmark input.txt
```

## Credits

- **Fiwa** — concept, design, architecture, every key insight
- **AI** — implementation
- **Yuta Mori** — [libdivsufsort](https://github.com/y-256/libdivsufsort) (MIT license), used for suffix array construction

## License

MIT
