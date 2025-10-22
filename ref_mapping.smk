# Load packages
import os                                  # Standard lib for paths/dirs if needed
import pandas as pd                        # Convenience; not strictly required in rules
import re                                  # 

# Global variables
BLAST = config["blast"]                    # BLAST parameter block
HITPOL = config["hit_policy"]              # Hit selection policy
WORK = config["work_dir"]                  # Output base directory
FROM = config["from_ref"]                  # Source reference key (e.g., "v8")
TO = config["to_ref"]                      # Target reference key (e.g., "v6")
FLANK = int(config["flank"])               # Flank size as integer
REFS = config["references"]                # Dict of references from config

# Prefer per-reference header rules; fall back to global keys if missing
ZPAD = int(REFS[FROM].get("header", {}).get("zpad", config.get("zpad", 0)))
CHROM_PREFIX = REFS[FROM].get("header", {}).get("prefix", config.get("chrom_prefix", ""))

# Regex to parse IDs into named groups chrom/pos
ID_REGEX = r"^(?:S|Chr|Chromosome|C|chr)?(?P<chrom>\d{1,3})_(?P<pos>\d+)$"

# Build a BLAST nucleotide database for the target reference
rule make_blastdb:
    input:
        fasta = lambda wc: REFS[TO]["fasta"]                # Target FASTA path
    output:
        f"{WORK}/blastdb/{TO}.nhr"                          # Presence of index file signals completion
    params:
        dbprefix = lambda wc: f"{WORK}/blastdb/{TO}"        # Output prefix for BLAST DB files
    conda:
        "envs/blast.yaml"                                   # Use the BLAST/py env
    shell:
        r"""
        mkdir -p {WORK}/blastdb                             # Ensure output dir exists.
        makeblastdb -in {input.fasta} -dbtype nucl -out {params.dbprefix}  # Create BLAST DB.
        """

# Parse SNP IDs (e.g., S07_2584157) into (chrom, pos), applying prefix/zpad
rule parse_ids:
    input:
        config["input_file"]                                # Raw SNP ID list
    output:
        tsv = f"{WORK}/parsed_ids.tsv"                      # Parsed table (ID, chrom, pos, status)
    params:
        regex = lambda wc: ID_REGEX,                        # Pass regex via lambda to avoid wildcard parsing
        chrom_prefix = lambda wc: CHROM_PREFIX,             # Source-specific chrom prefix
        zpad = lambda wc: ZPAD                              # Source-specific zero-padding
    conda:
        "envs/blast.yaml"                                   # Python + pandas environment
    script:
        "scripts/parse_ids.py"                              # Implements the parsing logic

# Extract flank sequences around each SNP from the source reference FASTA
rule extract_flanks:
    input:
        tsv = rules.parse_ids.output.tsv,                   # Parsed IDs with chrom/pos
        fasta = lambda wc: REFS[FROM]["fasta"]              # Source FASTA path
    output:
        fa = f"{WORK}/flanks/{FROM}.fa",                    # Fasta of all query windows
        map = f"{WORK}/flanks/{FROM}.map.tsv"               # Metadata: offsets and checks
    params:
        flank = FLANK                                       # Flank size
    conda:
        "envs/blast.yaml"                                   # Python + pyfaidx env
    script:
        "scripts/extract_flanks.py"                         # Extracts windows and writes FASTA/map

# BLAST the extracted flanks against the target reference DB
rule blast_flanks:
    input:
        fa = rules.extract_flanks.output.fa,                # Query sequences (source windows)
        db = rules.make_blastdb.output                      # Target BLAST DB (presence file)
    output:
        out = f"{WORK}/blast/{FROM}_to_{TO}.tsv"            # Tabular BLAST output (fmt 6)
    params:
        dbprefix = f"{WORK}/blastdb/{TO}",                  # Target DB prefix
        task = BLAST.get("task", "blastn"),                 # BLAST task to run
        evalue = BLAST["evalue"],                           # E-value threshold
        max_target_seqs = BLAST["max_target_seqs"],         # Max targets per query
        perc_identity = BLAST["perc_identity"],             # Identity threshold
        dust = BLAST.get("dust","no"),                      # Query masking
        threads = BLAST.get("threads", 4)                   # Threads for BLAST
    conda:
        "envs/blast.yaml"                                   # BLAST environment
    shell:
        r"""
        mkdir -p {WORK}/blast                               # Ensure output dir exists
        blastn -task {params.task} -query {input.fa} -db {params.dbprefix} \
               -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} \
               -perc_identity {params.perc_identity} -dust {params.dust} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
               -num_threads {params.threads} > {output.out} # Produce tabular results
        """

# Pick best hit per SNP and compute the lifted position on the target reference
rule pick_hit_and_convert:
    input:
        parsed = rules.parse_ids.output.tsv,                # Original parsed IDs
        fmap   = rules.extract_flanks.output.map,           # Window offsets to compute SNP position
        blast  = rules.blast_flanks.output.out              # BLAST alignments table
    output:
        csv = f"{WORK}/converted_coords.csv"                # Final joined table with statuses
    params:
        require_full_identity = HITPOL["full_identity"],    # Enforce 100% identity if true
        tie_break = HITPOL["tie_break"],                    # Strategy to resolve multiple hits
        to_ref = TO,                                        # Target reference label for output column
        from_ref = FROM                                     # Source reference label for output column
    conda:
        "envs/blast.yaml"                                   # Python + pandas environment
    script:
        "scripts/pick_hit.py"                               # Select best hit and compute coordinates