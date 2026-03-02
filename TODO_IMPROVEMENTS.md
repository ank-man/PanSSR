# PanSSR Improvement TODOs

## Phase 1 — Immediate fixes (implemented)
- [x] Avoid package-wide import failures by lazy-loading modules in `panssr/__init__.py`.
- [x] Improve genome/annotation pairing to use deterministic stem matching and avoid substring mismatches.
- [x] Reuse an open BAM handle during genotype runs to remove repeated open/close overhead.
- [x] Prevent genotype key collisions by indexing markers with `(chrom, start, end)` in memory/output.
- [x] Add `environment.yml` so CI workflow can resolve dependencies consistently.

## Phase 2 — Performance upgrades
- [ ] Parallelize genome-mode sequence processing with `ProcessPoolExecutor`.
- [ ] Precompile SSR regexes per motif length/threshold combination.
- [ ] Optimize ePCR candidate pairing using sorted reverse positions + binary search.
- [ ] Add optional chunked output writers for large runs.

## Phase 3 — Reliability/quality
- [ ] Add benchmark tests for SSR discovery and genotyping throughput.
- [ ] Add integration tests for genome/genotype modes using small fixtures.
- [ ] Add schema/version metadata to marker TSV and genotype CSV outputs.
- [ ] Improve error reporting with per-sample/per-locus failure summaries.

## Phase 4 — GUI roadmap
- [ ] Extract pipeline service layer from CLI (`run_genome_pipeline`, `run_genotype_pipeline`).
- [ ] Build FastAPI backend for job submission + status APIs.
- [ ] Build Streamlit GUI MVP (inputs, progress logs, table output, downloads).
- [ ] Add async task queue for long runs and persistent job history.
