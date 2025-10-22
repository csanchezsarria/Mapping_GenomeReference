# Load config
configfile: "config.yaml"

# Include modular rules for different stages of the workflow
include: "ref_mapping.smk"

# Final target: produce the lifted coordinates table.
rule all:
    input:
        f"{WORK}/converted_coords.csv"