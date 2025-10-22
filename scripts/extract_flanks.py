# Extract flanking windows around SNP coordinates from a source reference FASTA.
# Input:  parsed_ids.tsv (with columns: input_id, chrom, pos, status)
# Output: flanks FASTA (one record per SNP) + a map TSV with window metadata.

# Load packages
import pandas as pd                  # Tabular I/O and manipulation
from pyfaidx import Fasta            # Efficient random access to FASTA by name and slice
from pathlib import Path             # Convenient and safe path handling

# Snakemake-provided handles
tsv = snakemake.input.tsv            # Path to parsed_ids.tsv
fasta = snakemake.input.fasta        # Path to source reference FASTA (from_ref)
out_fa = snakemake.output.fa         # Path to write the flanking sequences FASTA
out_map = snakemake.output.map       # Path to write window metadata (TSV)
flank = int(snakemake.params.flank)  # Flank size on each side of the SNP (bp)

# Ensure output directory exists
Path(out_fa).parent.mkdir(parents = True, exist_ok = True)

# Load inputs
df = pd.read_csv(tsv, sep = "\t")    # Parsed IDs with columns from parse_ids.py
fa = Fasta(fasta, rebuild = False)   # Open the FASTA index; don't rebuild if .fai exists

records = []                         # FASTA records to write (strings with header+seq)
maps = []                            # Rows for the metadata TSV

# Iterate each parsed SNP
for i, r in df.iterrows():

    # Skip malformed IDs or missing chromosome names
    if r.get("status") != "ok" or pd.isna(r["chrom"]):
        continue

    chrom = str(r["chrom"])          # Chromosome name to query in FASTA
    pos = int(r["pos"])              # 1-based SNP coordinate

    # If chromosome not present in FASTA, mark and continue
    if chrom not in fa:
        maps.append(dict(input_id = r["input_id"], status = "no_chrom_in_src"))
        continue

    # Coordinates are 1-based in biological convention:
    # Extract [pos - flank, pos + flank], clipped to sequence bounds.
    start = max(1, pos - flank)
    end = min(len(fa[chrom]), pos + flank)

    # pyfaidx uses 0-based slices, end-exclusive; hence (start-1):end.
    seq = fa[chrom][start-1:end].seq.upper()

    # Build a descriptive query ID for traceability across steps
    qid = f"{r['input_id']}|{chrom}:{pos}|win:{start}-{end}"

    # Append a FASTA-formatted record
    records.append(f">{qid}\n{seq}\n")

    # Compute the SNP's offset within the window (0-based relative to 'start')
    snp_offset = pos - start

    # Record metadata for downstream coordinate conversion after BLAST
    maps.append(dict(
        input_id = r["input_id"],
        src_chrom = chrom,
        src_pos = pos,
        win_start = start,
        win_end = end,
        snp_offset = snp_offset,
        status = "ok"
    ))

# Write the flanking sequences FASTA
with open(out_fa, "w") as fh:
    fh.writelines(records)

# Write the window map/metadata TSV
pd.DataFrame(maps).to_csv(out_map, sep = "\t", index = False)
