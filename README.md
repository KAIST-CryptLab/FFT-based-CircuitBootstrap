# FFT-based CBS
This is an FFT-based implementation of ['FHEW-like Leveled Homomoprhic Evaluation: Refined Workflow and Polished Building Blocks'](https://eprint.iacr.org/2024/1318).

## Contents
We implement:
- benchmarks for
  - FFT-based circuit bootstrapping for the bitwise LHE (`bench_cbs`)
  - FFT-based circuit bootstrapping for the HP-LHE (`bench_hp_lhe`)
  - AES evaluation in LHE mode (`bench_aes`) and flex. LHE mode (`bench_aes_half_cbs`)

## How to use
- bench: `cargo bench --bench 'benchmark_name'`
  - Current sample size is set to 1000. It can be changed by modifying `config = Criterion::default().sample_size(1000);`
