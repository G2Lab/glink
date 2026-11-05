# Implementation of "GLINK: computing linkage disequilibrium from pangenome graphs"

## Overview

GLINK brings the familiar PLINK CLI to pangenome graphs. By matching PLINK's interface, GLINK can be easily swapped into existing pipelines whenever your input data shifts from VCFs to pangenome graphs.

Please email me at jb4816 [at] [columbia] [dot] edu if you have any questions!

## Build

```bash
cmake -B build && cmake --build build -j
```

## Dependencies

**Install these before building:**

- [CMake](https://cmake.org) 3.18 or newer
- C and C++17 compilers. GCC needs to be version 11 or newer.
- [vg](https://github.com/vgteam/vg) on the `PATH`. Used here with v1.73.0 and v1.78.0.
- [zlib](https://zlib.net), with headers, to read GFAs and hash them
- [OpenSSL](https://www.openssl.org), with headers, for SHA-256

**Fetched or included automatically:**

- [HDF5](https://github.com/HDFGroup/hdf5) 1.14.6 and [HighFive](https://github.com/highfive-devs/highfive) 3.1.1
  are downloaded by `cmake` and linked statically.
- [nlohmann/json](https://github.com/nlohmann/json) 3.12.0 and [GoogleTest](https://github.com/google/googletest)
  are in the repository.

## Use

```bash
build/bin/glink --gfa chr1.gfa --out-prefix idx/chr1 --filter @nodes.txt --matrix-type dosage \
          --condition-on-parents --collapse-backbone-parents --filter common --no-collapse \
          --samples cohort.txt --threads 8 \
          --matrix-out m --r2-out r2.ld --dprime-out dprime.ld --blocks
```

r² and D′ are written in PLINK 1.9's exact formats (see [File formats](#file-formats)):

| Option | File | Same as |
| --- | --- | --- |
| `--r2-out r2.ld` | `r2.ld`: the r² matrix, tab-separated; `r2.snplist`: its row IDs in order | `plink --r2 square --write-snplist` |
| `--dprime-out d.ld` | `d.ld`: one row per pair, `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2 DP` | `plink --r2 inter-chr dprime` |

A name ending in `.ld.bin` writes raw little-endian float32 instead, `n` x `n`, row by row, with the same
`.snplist`: much faster than text, and it keeps every digit of r². With `--dprime-out` it holds D′ alone.

D′ is signed for phased or presence matrices (as PLINK's `dprime-signed`) and a genotype-likelihood
\|D′\| for dosage. PLINK has no matrix form of D′, so its table is the only form.

`--threads` also parallelises r² and D′ (OpenMP).

## Standalone runs and `glink index`

`--out-prefix` makes one command do everything from a GFA. It builds the index on the first
run and reuses it on later ones:

```bash
build/bin/glink --gfa chr1.gfa --out-prefix idx/chr1 --condition-on-parents --filter @nodes.txt \
          --threads 8 --r2-out r2.ld --blocks --blocks-out blocks.tsv
```

`glink index` builds the same files without the analysis:

```bash
build/bin/glink index --gfa chr1.gfa --out-prefix idx/chr1 --threads 8
```

| File | Contents |
| --- | --- |
| `PREFIX.glink/graph.h5` | The parsed GFA: node x haplotype presence bits, links, and the reference and backbone walks |
| `PREFIX.glink/tree.h5` | Conditioning tree: `vg snarls` rooted on `--backbone`, each node's innermost snarl, and so each node's parent (that snarl's start node) |
| `PREFIX.glink/coords.h5` | Node coordinates on `--coordinate-reference` |

A prefix ending in `/` names the index after the GFA.

* **Resuming:** each part is checked on its own, so an interrupted build redoes only what is
  missing. Only building the tree needs `vg`. A run that reuses an index does not.
* **Without `--out-prefix`**, the parsed graph is saved under `.glink-cache/graph/` and loaded instead of
  parsing the GFA again, until the GFA's contents change.
* **`--rebuild`** replaces the whole index.
* **Stale files are errors.** A file that does not match the GFA's contents (SHA-256 of the
  uncompressed text) or the backbone stops the run and names `--rebuild`. Nothing is replaced unasked.
  Every run hashes the whole GFA to check this, and `--gfa` stays required. Moving, copying or
  touching the GFA or an index does not invalidate anything.
  A different `vg` version does not invalidate an existing tree.
* **`--split-inversions`** keeps its own `graph_oriented.h5`, because `graph.h5` has no orientation. On an index
  that has only `graph.h5`, the first such run reads the GFA once more.
* `--out-prefix` replaces `--node-coordinates`, `--conditioning-tree`, `--presence-bundle` and the other
  cache options, and cannot be combined with them. Without it, `glink` behaves as before.

## Tests

```bash
cmake --build build -j --target glink-tests && build/bin/glink-tests
```