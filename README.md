# Ref mapping by BLAST with Snakemake 🧬

Snakemake pipeline to *lift/map* SNP coordinates from a **source** reference to a **target** reference using flanking windows and **BLAST+**. Given SNP IDs like `C01_7120920`, the workflow extracts a flanking window around each SNP from the source FASTA, runs **blastn** against the target reference, and selects the best hit to report the lifted coordinate on the target genome. The step‑by‑step logic is implemented in the scripts under `scripts/`. 




## Authors 🙋

For questions or feedback about this pipeline, please contact:

- **Camilo E. Sánchez-Sarria** - [CGIAR]: c.e.sanchez@cgiar.org


## Repository layout

.
├── config.yaml
├── envs
│ └── blast.yaml
├── input_snps.txt
├── ref_mapping.smk
├── results/
├── scripts
│ ├── extract_flanks.py
│ ├── parse_ids.py
│ └── pick_hit.py
└── Snakefile


## Workflow overview 🚀

input_snps.txt
└─ parse_ids.py ───────────→ parsed_ids.tsv
└─ extract_flanks.py (+ from_ref.fasta) → flanks.fa + window_map.tsv
└─ BLASTN (flanks.fa vs to_ref.fasta) → blast.tsv (fmt 6)
└─ pick_hit.py → results/mapping.csv

- **ID parsing** — takes `input_snps.txt`, applies a configurable regex with named groups `(?P<chrom>...)` and `(?P<pos>...)`, builds chromosome names using a prefix and zero‑padding if needed, and writes `parsed_ids.tsv`. Unparseable IDs are flagged with `status=bad_id`. :contentReference[oaicite:1]{index=1}  
- **Flank extraction** — for each valid ID, opens the **source** FASTA with `pyfaidx`, extracts a window `[pos - flank, pos + flank]` (clipped to contig bounds), emits a query FASTA and a `window_map.tsv` with metadata (including `snp_offset`). Chromosomes missing from the source FASTA are flagged as `no_chrom_in_src`. :contentReference[oaicite:2]{index=2}  
- **BLASTn** — runs **BLAST+** (outfmt 6) of the flanks against the **target** reference, producing `blast.tsv` with standard columns.  
- **Hit selection / lift** — filters candidates (optionally requiring 100% identity), applies a tiebreak policy, computes `output_pos` using `snp_offset`, determines `strand`, and writes `results/mapping.csv` with:
  `input_id, from_ref, to_ref, output_id, output_chrom, output_pos, strand, identity, aln_len, evalue, bitscore, status`. :contentReference[oaicite:3]{index=3}



## Installation and dependencies 🛠️




## Usage ▶️




## Citing this pipeline 📌




## License 📄