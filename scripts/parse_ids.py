# Parse SNP IDs into (chromosome, position) using a configurable regex.
# Input : a 1-column file with IDs (CSV or TSV), e.g., "S07_2584157"
# Output: parsed_ids.tsv with columns: input_id, chrom, pos, status

# Load packages
import re, sys, pandas as pd                    # regex for parsing, pandas for tabular I/O

# Snakemake-provided handles
inp = snakemake.input[0]                        # Path to the raw ID list
out = snakemake.output.tsv                      # Path to write the parsed table (TSV)
regex = re.compile(snakemake.params.regex)      # Compiled regex with named groups (?P<chrom>, ?P<pos>)
chrom_prefix = snakemake.params.chrom_prefix    # Prefix to build FASTA-style chrom names
zpad = int(snakemake.params.zpad)               # Zero-padding width for chromosome numbers

# Read a single-column file of IDs; try CSV first, then fallback to TSV/space-delimited
try:
    df = pd.read_csv(inp, header = None)
except:
    df = pd.read_csv(inp, header=None, sep="\t")

df.columns = ["input_id"]                       # Normalize the column name
rows = []                                       # Accumulator for parsed rows

for idv in df["input_id"].astype(str):
    s = idv.strip()                             # Trim whitespace/newlines
    m = regex.match(s)                          # Full-string match (use ^...$ in regex for safety)
    if not m:

        # Could not parse: record as a bad ID
        rows.append(dict(input_id = s, chrom = None, pos = None, status = "bad_id"))
        continue

    chrom_num = m.group("chrom")                # Extract chromosome number as text
    pos = int(m.group("pos"))                   # Extract position as integer

    # Optional zero-padding of chromosome number (e.g., "7" -> "07")
    if zpad > 0:
        chrom_num = str(chrom_num).zfill(zpad)

    # Build the final chromosome name expected by the source FASTA (prefix may be empty)
    chrom_name = f"{chrom_prefix}{chrom_num}" if chrom_prefix else str(chrom_num)

    # Record a successful parse
    rows.append(dict(input_id = s, chrom = chrom_name, pos = pos, status = "ok"))

# Write the parsed table as TSV for downstream steps
pd.DataFrame(rows).to_csv(out, sep = "\t", index=False)