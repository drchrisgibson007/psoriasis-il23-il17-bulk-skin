# Measurement limitations of IL-23/IL-17 bulk-skin stratification in psoriasis

Analysis code for a manuscript targeting the **Journal of Dermatological Science**. The scripts score lesional bulk skin on an IL-23-axis module and an IL-17-response module across six public cohorts, then repeat that stratification on deposited values, on RMA re-normalized arrays, and after QC removal.

**Authors:** Christopher G. Gibson, DO; Alyssa Jack, PA-C; Jared Gibson; Liam Smith; Illari Research Institute

## How to run

From the repository root:

```bash
bash run_all.sh
```

`run_all.sh` sources `config.sh` (sets `GRG_ROOT`) and runs `analysis/run_all.sh`. That driver performs three passes:

1. **Deposited values.** `analysis/scripts/analyze.py` and `analysis/scripts/sens.py` read submitter-processed microarray matrices and deposited RNA-seq counts or RPKM, and write `analysis/out/`.
2. **RMA re-normalization.** `analysis/scripts/analyze_rma.py` and `analysis/scripts/sens_rma.py` repeat the stratification on RMA matrices and write `analysis/out_rma/`. CEL files are summarized with `analysis/scripts/rma.R`; batch diagnostics are in `analysis/scripts/batch.R`.
3. **QC removal.** `QCDROP=1 python3 analyze_rma.py` drops arrays flagged by `analysis/scripts/qcflags.py` (NUSE and RLE) and writes `analysis/out_rma_qc/`. This pass stops after those outputs: three QC-dropped GSE117468 baseline arrays are absent from the response section, which this sensitivity analysis does not use.

`analysis/scripts/overlap_check.py` records cohort overlap from GEO sample metadata.

Python 3 with numpy, pandas, scipy, and matplotlib is required for `bash run_all.sh`. Re-normalization from CEL files also needs R with the Bioconductor packages affy, affyPLM, limma, and sva.

Place series matrices and RNA-seq supplements in `analysis/data/`, and RMA matrices plus QC tables in `analysis/rma/`. The comparison at the end of `analysis/run_all.sh` reads reference tables from `analysis/outputs_reference/`. Those files are public GEO inputs and local results. They are not part of this repository. Fetch commands are below.

`python3 make_figures.py` redraws Figures 1 and 2 from a lesion-level scores table. Point `SRC` at the top of that script to the `lesion_level_scores.tsv` produced by a run.

## Public GEO data

Expression inputs are public Gene Expression Omnibus series. This repository holds the analysis scripts, the CEL checksum manifest, and the series-matrix MD5 list. Large CEL archives and expression matrices are not included.

| Accession | Role |
| --- | --- |
| GSE13355 | GPL570 microarray |
| GSE14905 | GPL570 microarray |
| GSE30999 | GPL570 microarray |
| GSE117468 | GPL570 microarray |
| GSE121212 | RNA-seq read counts |
| GSE54456 | RNA-seq RPKM |

Series matrices and supplementary count tables (about 300 MB):

```bash
bash analysis/raw_rebuild/fetch_full_inputs.sh
```

CEL archives for the four microarray cohorts (about 3.6 GB), checked against MD5 sums in `analysis/raw_rebuild/raw_manifest.tsv`:

```bash
bash analysis/raw_rebuild/fetch_raw.sh /path/to/dest
```

MD5 checksums for the series-matrix bundle, as downloaded on 2026-10-08, are in `analysis/raw_rebuild/full_inputs_md5_as_downloaded_2026-10-08.txt`. GEO FTP URLs are in the manifest and in `analysis/raw_rebuild/fetch_full_inputs.sh`.

## License

MIT. See [LICENSE](LICENSE).
