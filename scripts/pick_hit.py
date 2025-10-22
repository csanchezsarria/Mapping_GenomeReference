# Select the best BLAST hit per SNP and compute the lifted coordinate on the target reference.
# Final columns (as requested):
#   input_id, from_ref, to_ref, output_id, output_chrom, output_pos, strand,
#   identity, aln_len, evalue, bitscore, status

# Load packages
import pandas as pd                # Tabular I/O
import numpy as np                 # NaN handling and vector ops

# Snakemake-wired inputs
parsed = pd.read_csv(snakemake.input.parsed, sep = "\t")  # From parse_ids.py
fmap   = pd.read_csv(snakemake.input.fmap,   sep = "\t")  # From extract_flanks.py
blast  = pd.read_csv(snakemake.input.blast, # Tabular BLAST fmt 6
    sep = "\t",
    header = None,
    names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
    "qend", "sstart", "send", "evalue", "bitscore"],
)

# Snakemake-wired params
req_full = bool(snakemake.params.require_full_identity)  # If True, enforce 100% identity
tie_break = snakemake.params.tie_break                   # Hit selection strategy
to_ref = snakemake.params.to_ref                         # Target reference label (for output)
from_ref = snakemake.params.from_ref                     # Source reference label (for output)

# Keep only the columns we need from 'parsed' and 'fmap' (and give unique names for statuses)
parsed_min = parsed[["input_id","status"]].rename(columns = {"status": "parsed_status"})
fmap_min = fmap[["input_id","snp_offset","status"]].rename(columns = {"status": "fmap_status"})

# qseqid layout is "input_id|chrom:pos|win:start-end"; extract the leading input_id token
blast["input_id"] = blast["qseqid"].str.split("|", n = 1, expand = True)[0]

# Optionally require perfect identity (your "Phytozome-like" strict mode)
if req_full:
    blast = blast[blast["pident"] == 100.0].copy()

# Rank candidate hits per input_id according to the requested tiebreak policy
if tie_break == "bitscore_then_len":
    # Prefer highest bitscore, then longest alignment
    blast = blast.sort_values(["input_id","bitscore","length"], ascending = [True, False, False])
elif tie_break == "longest":
    # Prefer longest alignment, then bitscore
    blast = blast.sort_values(["input_id","length","bitscore"], ascending = [True, False, False])
elif tie_break == "leftmost":
    # Prefer leftmost target start (min(sstart, send)), then bitscore
    left = blast.copy()
    left["ss"] = left[["sstart","send"]].min(axis = 1)
    blast = left.sort_values(["input_id","ss","bitscore"], ascending = [True, True, False])
else:
    # Default: highest bitscore wins
    blast = blast.sort_values(["input_id","bitscore"], ascending = [True, False])

# Take the top-ranked hit per input_id (if any)
best = blast.groupby("input_id").head(1).copy()

# Attach SNP offset and source QC status (no need to bring src coordinates)
best = best.merge(fmap_min, on = "input_id", how="left")

# Convert SNP offset within the window to a target genomic coordinate
# For + strand: target_pos = min(sstart, send) + snp_offset
# For - strand: target_pos = max(sstart, send) - snp_offset
def compute_target_pos(row):
    if pd.isna(row.get("snp_offset")):
        return np.nan

    off = int(row["snp_offset"])

    if row["sstart"] <= row["send"]:  # '+' strand
        return int(min(row["sstart"], row["send"]) + off)
    else:                              # '-' strand
        return int(max(row["sstart"], row["send"]) - off)

# Derive strand, target (output) position and chromosome from the chosen hit
best["strand"] = np.where(best["sstart"] <= best["send"], "+", "-")
best["output_pos"] = best.apply(compute_target_pos, axis = 1)
best["output_chrom"] = best["sseqid"]

# Keep/rename the hit metrics we want to expose
hits = best[["input_id", "output_chrom", "output_pos", "strand", "pident", "length", 
    "evalue", "bitscore", "fmap_status"]].copy()
hits.rename(columns = {"pident": "identity", "length": "aln_len"}, inplace = True)

# Build "output_id" as "<output_chrom>_<output_pos>"; leave NaN if no position
hits["output_id"] = np.where(
    hits["output_pos"].notna(),
    hits["output_chrom"].astype(str) + "_" + hits["output_pos"].astype("Int64").astype(str),
    np.nan
)

# Merge back to keep all original IDs (including bad_id / no_chrom_in_src / no_hit)
out = parsed_min.merge(hits, on = "input_id", how = "left")

# Direction labels
out["from_ref"] = from_ref
out["to_ref"] = to_ref

# Status synthesis (pipeline-level QC):
# ok → we computed output_pos
# bad_id → regex failed in parse_ids
# no_chrom_in_src → source chrom missing in source FASTA (prefix/padding issue)
# no_hit → no BLAST hit survived filters / tiebreak
out["status"] = np.where(
    out["output_pos"].notna(), "ok",
    np.where(
        out["parsed_status"].eq("bad_id"), "bad_id",
        np.where(out["fmap_status"].eq("no_chrom_in_src"), "no_chrom_in_src", "no_hit")
    )
)

# Final column order (no src_chrom/src_pos)
final_cols = ["input_id", "from_ref", "to_ref", "output_id", "output_chrom", "output_pos", "strand",
    "identity", "aln_len", "evalue", "bitscore", "status"
]

# Write final CSV
out[final_cols].to_csv(snakemake.output.csv, index = False)