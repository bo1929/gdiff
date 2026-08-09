# Tests

## Unit tests

Source files live in `unit/`. Built with [doctest](https://github.com/doctest/doctest) (header in `vendor/doctest/`).

```
make test            # build & run all tests
```

The test binary `build/test_gdiff` links against all project objects except `gdiff.o` (main).

## Regression tests

```
make test-unit       # unit tests
make test-regression # sketch + map on all genome pairs, diff against gt/
make test            # unit + regression
```

Ground truth lives in `gt/` with two suffixes:
- `*.enum.txt` — `--enum-only` mode output
- `*.cont.txt` — continuous (default) mode output

On the first run without `*.cont.txt` files, the script generates them automatically.

## Test data

- `genomes/` — compressed FASTA files (`.fna.gz`)
- `genome_names.txt` — one genome ID per line
- `genome_pairs.txt` — tab-separated query/ref pairs
