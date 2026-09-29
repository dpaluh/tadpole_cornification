"""
Build type I and type II keratin FASTAs for phylogenetic analysis.
Headers: >Species_code|gene|protein_id
"""

from pathlib import Path
import pandas as pd

# Inputs
KERATINS_TSV = "/project/inbreh/frog_keratin/hmm_processed/keratin_genes_classified.tsv"
NCBI_TSV = "/project/inbreh/frog_keratin/subs_ref_genomes/ncbi_dataset.tsv"
PROTEOME_DIR = Path("/project/inbreh/frog_keratin/all_ref_genomes/protein_fas")
OUT_DIR = Path("/project/inbreh/frog_keratin/tree_input")
OUT_DIR.mkdir(exist_ok=True)


def make_species_code(organism_name):
    """
    Turn 'Rana temporaria' into 'Ran_tem', 'Gastrophryne carolinensis' into 'Gas_car'.
    First 3 letters of genus + first 3 of species.
    """
    parts = organism_name.split()
    if len(parts) < 2:
        return organism_name.replace(" ", "_")
    return parts[0][:3].capitalize() + "_" + parts[1][:3].lower()


# --------------------------------------------------------------------------
# 1. Build assembly accession -> species code mapping
# --------------------------------------------------------------------------
ncbi = pd.read_csv(NCBI_TSV, sep="\t")

# Some rows have both a GCF and its paired GCA; we want to map both to the same code
accession_to_code = {}
accession_to_organism = {}
for _, row in ncbi.iterrows():
    org = row["Organism Name"]
    if pd.isna(org):
        continue
    code = make_species_code(org)
    for col in ["Assembly Accession", "Assembly Paired Assembly Accession"]:
        acc = row.get(col)
        if pd.notna(acc):
            accession_to_code[acc] = code
            accession_to_organism[acc] = org

# Add mouse and Anolis manually since they aren't in ncbi_dataset.tsv
manual_additions = {
    "GCF_000001635.27": ("Mus_mus", "Mus musculus"),
    "GCF_035594765.1":  ("Ano_car", "Anolis carolinensis"),
}
for acc, (code, org) in manual_additions.items():
    accession_to_code.setdefault(acc, code)
    accession_to_organism.setdefault(acc, org)

print(f"Species mapping built: {len(set(accession_to_code.values()))} unique species codes")


# --------------------------------------------------------------------------
# 2. Load keratins and add species codes
# --------------------------------------------------------------------------
keratins = pd.read_csv(KERATINS_TSV, sep="\t")

def species_to_code(species_str):
    # species column looks like: "GCF_000001635.27_GRCm39"
    # extract the accession (first two underscore-separated parts)
    acc = "_".join(species_str.split("_")[:2])
    return accession_to_code.get(acc, acc)

keratins["species_code"] = keratins["species"].apply(species_to_code)

# Check we mapped everything
unmapped = keratins[keratins["species_code"].str.startswith("GC")]
if len(unmapped) > 0:
    print(f"WARNING: {len(unmapped)} keratins from unmapped species:")
    print(unmapped["species"].unique())

print("Keratin counts per species code:")
print(keratins.groupby(["species_code", "classification"]).size().unstack(fill_value=0))


# --------------------------------------------------------------------------
# 3. Build a species -> proteome file mapping
# --------------------------------------------------------------------------
species_to_proteome = {}
for faa in PROTEOME_DIR.glob("*.faa"):
    # File names like "GCF_000001635.27_GRCm39_protein.faa"
    accession = "_".join(faa.name.split("_")[:2])
    species_to_proteome[accession] = faa

# Map from the dataframe's species column to proteome files
species_col_to_proteome = {}
for sp in keratins["species"].unique():
    acc = "_".join(sp.split("_")[:2])
    if acc in species_to_proteome:
        species_col_to_proteome[sp] = species_to_proteome[acc]
    else:
        print(f"WARNING: no proteome for {sp}")

print(f"Proteomes mapped: {len(species_col_to_proteome)}/{keratins['species'].nunique()}")


# --------------------------------------------------------------------------
# 4. Extract sequences and write FASTAs
# --------------------------------------------------------------------------
def read_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file."""
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, "".join(seq)


def extract_protein_id(header):
    """Get the protein ID from a FASTA header. First whitespace-delimited token."""
    return header.split()[0]


# Group keratins by proteome file so we only read each proteome once
records_to_extract = {}  # proteome_path -> {protein_id: (species_code, gene, classification)}
for _, row in keratins.iterrows():
    proteome = species_col_to_proteome.get(row["species"])
    if proteome is None:
        continue
    records_to_extract.setdefault(proteome, {})[row["target"]] = (
        row["species_code"],
        row["gene"],
        row["classification"],
    )

# Now scan each proteome and pull the sequences
extracted = {"typeI": [], "typeII": []}  # list of (header, seq) per type
missing = []

for proteome, wanted in records_to_extract.items():
    print(f"Scanning {proteome.name}...", flush=True)
    found_ids = set()
    for header, seq in read_fasta(proteome):
        pid = extract_protein_id(header)
        if pid in wanted:
            species_code, gene, classification = wanted[pid]
            new_header = f"{species_code}|{gene}|{pid}"
            extracted[classification].append((new_header, seq))
            found_ids.add(pid)
    # Log any wanted IDs we didn't find
    for pid in set(wanted.keys()) - found_ids:
        missing.append((proteome.name, pid))

if missing:
    print(f"\nWARNING: {len(missing)} keratins not found in proteome files:")
    for prot, pid in missing[:10]:
        print(f"  {prot}: {pid}")


# --------------------------------------------------------------------------
# 5. Write outputs
# --------------------------------------------------------------------------
for classification, records in extracted.items():
    out_path = OUT_DIR / f"keratins_{classification}.fasta"
    with open(out_path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            # Wrap sequence at 60 chars per line
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    print(f"Wrote {len(records)} sequences to {out_path}")