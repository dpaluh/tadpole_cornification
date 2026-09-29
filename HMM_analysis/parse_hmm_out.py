# imports
from pathlib import Path
import pandas as pd
import re
import subprocess
import gzip



# Set directory to hmm searches
TBL_DIR = Path("/project/inbreh/frog_keratin/hmm_searches")

# parse all of the output
def parse_tblout(path):
    """Parse hmmsearch --tblout output into a list of dicts."""
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 18)
            if len(parts) < 19:
                continue
            hits.append({
                "target": parts[0],
                "evalue": float(parts[4]),
                "score": float(parts[5]),
                "description": parts[18].rstrip(),
            })
    return hits

# Load everything into one dataframe
records = []
for tbl in sorted(TBL_DIR.glob("*.tbl")):
    name = tbl.stem  # e.g. GCA_027917425.1_aGasCar1.pri_genomic_typeI
    if name.endswith("_typeI"):
        species = name[:-len("_typeI")]
        hmm_type = "typeI"
    elif name.endswith("_typeII"):
        species = name[:-len("_typeII")]
        hmm_type = "typeII"
    else:
        continue
    for hit in parse_tblout(tbl):
        hit["species"] = species
        hit["hmm"] = hmm_type
        records.append(hit)

df = pd.DataFrame(records)
print(df.shape)
print(df.head())
print(df["species"].value_counts())


# Pivot to one row per (species, target) with both HMM scores
wide = df.pivot_table(
    index=["species", "target", "description"],
    columns="hmm",
    values="score",
    aggfunc="max",  # in case of duplicates, keep highest
).reset_index()

# Fill missing scores with 0 (protein hit by only one HMM)
wide["typeI"] = wide["typeI"].fillna(0)
wide["typeII"] = wide["typeII"].fillna(0)

# Also pull in the best e-value for each protein (lowest = best)
evalue_min = df.groupby(["species", "target"])["evalue"].min().reset_index()
evalue_min.columns = ["species", "target", "best_evalue"]
wide = wide.merge(evalue_min, on=["species", "target"])

# Classify by reciprocal best score
def classify(row):
    if row["typeI"] > row["typeII"]:
        return "typeI"
    elif row["typeII"] > row["typeI"]:
        return "typeII"
    else:
        return "tie"

wide["classification"] = wide.apply(classify, axis=1)
wide["score_diff"] = (wide["typeI"] - wide["typeII"]).abs()

print(wide.shape)
print(wide.head())
print(wide["classification"].value_counts())



# Remove the non-keratin genes:

EXCLUDE_PATTERNS = [
    "vimentin", "desmin", "glial fibrillary", "alpha-internexin",
    "internexin", "peripherin", "neurofilament", "lamin", "prelamin",
    "nestin", "synemin", "syncoilin", "phakinin", "filensin",
    "intermediate filament family orphan",
    "intermediate filament protein ON3",
    "low molecular weight neuronal",
    "non-homologous end joining factor IFFO1",
    "thread biopolymer filament",
    "ouroboros", "beaded filament structural",
    "lamin tail domain",
]
exclude_re = re.compile("|".join(EXCLUDE_PATTERNS), re.IGNORECASE)

wide["is_other_IF"] = wide["description"].str.contains(exclude_re, regex=True)
wide["is_uncharacterized"] = wide["description"].str.lower().str.startswith("uncharacterized protein") \
                             | wide["description"].str.lower().str.startswith("low quality protein: uncharacterized")

# Summary
print("Total:", len(wide))
print("Other IFs (excluded):", wide["is_other_IF"].sum())
print("Uncharacterized:", wide["is_uncharacterized"].sum())
print("Likely keratins:", (~wide["is_other_IF"] & ~wide["is_uncharacterized"]).sum())

# Counts per species after filtering
keratins = wide[~wide["is_other_IF"] & ~wide["is_uncharacterized"]]
print("\nKeratins per species/type:")
print(keratins.groupby(["species", "classification"]).size().unstack(fill_value=0))


# Distribution of score differences for the classified hits
print("Score difference distribution for classified keratins:")
print(keratins["score_diff"].describe())
print("\nKeratins with small score margin (potentially ambiguous):")
ambiguous = keratins[keratins["score_diff"] < 50]
print(f"N = {len(ambiguous)}")
print(ambiguous[["species", "target", "typeI", "typeII", "classification", "description"]].head(20))


##############################################################################
#    Get gene to protein mapping from GFF file to collapise isoforms
##############################################################################

GFF_DIR = Path("/project/inbreh/frog_keratin/all_ref_genomes/gffs")

def parse_gff_for_proteins(gff_path):
    records = {}
    open_fn = gzip.open if str(gff_path).endswith(".gz") else open
    with open_fn(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            start, end = int(start), int(end)
            attr_dict = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v
            pid = attr_dict.get("protein_id")
            if not pid:
                continue
            # Try gene, then locus_tag (GCA GFFs), then Name, then fall back to protein ID
            gene = attr_dict.get("gene") or attr_dict.get("locus_tag") or attr_dict.get("Name", pid)
            if pid not in records:
                records[pid] = {
                    "gene": gene,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                }
            else:
                records[pid]["start"] = min(records[pid]["start"], start)
                records[pid]["end"] = max(records[pid]["end"], end)
    return records



# Map species ID from dataframe to GFF file
species_to_gff = {}
for gff in GFF_DIR.glob("*.gff*"):
    # Use just the accession (first two underscore-separated parts) as the key
    accession = "_".join(gff.name.split("_")[:2])
    species_to_gff[accession] = gff

# Build a mapping from df species name to GFF using accession
df_to_gff = {}
for sp in wide["species"].unique():
    accession = "_".join(sp.split("_")[:2])
    if accession in species_to_gff:
        df_to_gff[sp] = species_to_gff[accession]
    else:
        print(f"No GFF for {sp}")

print(f"Matched {len(df_to_gff)}/{len(wide['species'].unique())} species")

# Parse all GFFs and build a comprehensive protein-to-location mapping
all_records = []
for sp, gff_path in df_to_gff.items():
    print(f"Parsing {sp}...", flush=True)
    records = parse_gff_for_proteins(gff_path)
    for pid, info in records.items():
        all_records.append({
            "species": sp,
            "target": pid,
            **info,
        })
    print(f"  {len(records)} proteins")

protein_locations = pd.DataFrame(all_records)
print(f"\nTotal protein records: {len(protein_locations)}")
print(protein_locations.head())






# Merge location info into the keratin hits
keratins_loc = keratins.merge(
    protein_locations[["species", "target", "gene", "chrom", "start", "end", "strand"]],
    on=["species", "target"],
    how="left",
)

# Check whether all keratin hits got location info
missing_loc = keratins_loc[keratins_loc["gene"].isna()]
print(f"Keratin hits without location info: {len(missing_loc)}")
if len(missing_loc) > 0:
    print(missing_loc[["species", "target", "description"]].head(10))

# Collapse to one row per gene, keeping the best-scoring isoform
# Use the max of (typeI, typeII) score as the score to compare on
keratins_loc["best_score"] = keratins_loc[["typeI", "typeII"]].max(axis=1)
genes = (
    keratins_loc.sort_values("best_score", ascending=False)
    .drop_duplicates(subset=["species", "gene"], keep="first")
    .reset_index(drop=True)
)

print(f"\nBefore collapsing: {len(keratins_loc)} protein-level hits")
print(f"After collapsing:  {len(genes)} gene-level hits")

print("\nGene-level keratin counts per species/type:")
genes_per_spec_parsed = genes.groupby(["species", "classification"]).size().unstack(fill_value=0)
print(genes_per_spec_parsed)

genes_per_spec_parsed.to_csv("/project/inbreh/frog_keratin/hmm_processed/keratin_genes_per_species_classified.tsv", sep="\t", index=True)

genes.to_csv("/project/inbreh/frog_keratin/hmm_processed/keratin_genes_classified.tsv", sep="\t", index=False)
print(f"Saved {len(genes)} gene-level keratin hits")



# For each species/classification, show which chromosomes the keratins fall on
cluster_check = (
    genes.groupby(["species", "classification", "chrom"])
    .size()
    .reset_index(name="n_genes")
    .sort_values(["species", "classification", "n_genes"], ascending=[True, True, False])
)

# Print one species per group to see the distribution
print("Chromosome distribution of keratin hits per species/type:")
for sp in sorted(genes["species"].unique()):
    print(f"\n{sp}")
    sub = cluster_check[cluster_check["species"] == sp]
    for _, row in sub.iterrows():
        print(f"  {row['classification']:7s} {row['chrom']:25s} n={row['n_genes']}")


# save it to tsv:
cluster_check.to_csv("/project/inbreh/frog_keratin/hmm_processed/chrom_distribution_per_species.tsv", sep="\t", index=False)




# Most keratin hits cluster on 2 chromosomes in each genome
#    in each, there is one Type I that is on the Type II chromosome, should be Krt 18
#     in most, there are a few other scattered hits, need to check on those, too
#     X. laevis has type I and type II each in large numbers on 2 chromosomes, as expected because it is tetraploid


# Exclude X. laevis from the scattered analysis because it is tetraploid
non_laevis = genes[genes["species"] != "GCF_017654675.1_Xenopus_laevis_v10.1_genomic"]

# Get the pimrary clusters - i.e., where most type I are and where most type II are
primary_clusters = (
    non_laevis.groupby(["species", "classification", "chrom"])
    .size()
    .reset_index(name="n_genes")
    .sort_values(["species", "classification", "n_genes"], ascending=[True, True, False])
    .drop_duplicates(subset=["species", "classification"], keep="first")
    [["species", "classification", "chrom"]]
    .rename(columns={"chrom": "primary_chrom"})
)

genes_flagged = non_laevis.merge(primary_clusters, on=["species", "classification"], how="left")
genes_flagged["in_primary_cluster"] = genes_flagged["chrom"] == genes_flagged["primary_chrom"]

# get the ones that are scattered to other places
scattered = genes_flagged[~genes_flagged["in_primary_cluster"]].sort_values(
    ["species", "classification", "chrom"]
)
print(f"Scattered (non-primary cluster) hits: {len(scattered)}")
print(scattered[["species", "classification", "chrom", "gene", "description"]].to_string())

# Save scatterd hits to look at later
scattered[["species", "classification", "chrom", "gene", "description",
           "typeI", "typeII", "start", "end", "strand", "best_evalue"]].to_csv( "/project/inbreh/frog_keratin/hmm_processed/scattered_hits.tsv", sep="\t", index=False
)



### Check for anchor genes that should flank each of the keratin clusters

ANCHOR_GENES = {
    "typeI": ["SMARCE1", "EIF1"],
    "typeII": ["FAIM2", "EIF4B"],
}

def find_anchors_in_gff(gff_path, anchor_names):
    """
    Find anchor genes in a GFF by matching gene names (case-insensitive).
    Returns dict: {anchor_name: [(chrom, start, end, strand, full_gene_name), ...]}
    """
    # Build case-insensitive regex for anchor matching on the gene= attribute
    open_fn = gzip.open if str(gff_path).endswith(".gz") else open
    found = {a: [] for a in anchor_names}
    upper_targets = {a.upper(): a for a in anchor_names}
    
    with open_fn(gff_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            attr_dict = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v
            gene = attr_dict.get("gene") or attr_dict.get("Name", "")
            if gene.upper() in upper_targets:
                key = upper_targets[gene.upper()]
                found[key].append((chrom, int(start), int(end), strand, gene))
    return found

# Find anchors in all species
all_anchors = []
for sp, gff_path in df_to_gff.items():
    print(f"Searching anchors in {sp}...", flush=True)
    anchors = find_anchors_in_gff(gff_path, ["SMARCE1", "EIF1", "FAIM2", "EIF4B"])
    for anchor, hits in anchors.items():
        for chrom, start, end, strand, gene_name in hits:
            all_anchors.append({
                "species": sp,
                "anchor": anchor,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_name": gene_name,
            })

anchors_df = pd.DataFrame(all_anchors)
print(f"\nTotal anchor hits: {len(anchors_df)}")
print(anchors_df["anchor"].value_counts())


## missing anchors in several species, let's see what's going on:

# Per-species summary of which anchors were found
species_anchor_matrix = (
    anchors_df.groupby(["species", "anchor"]).size().unstack(fill_value=0)
)
print("Number of anchor hits per species:")
print(species_anchor_matrix)

# Identify species missing each anchor
all_species = set(df_to_gff.keys())
for anchor in ["SMARCE1", "EIF1", "FAIM2", "EIF4B"]:
    found_in = set(anchors_df[anchors_df["anchor"] == anchor]["species"])
    missing = all_species - found_in
    print(f"\n{anchor} missing in {len(missing)} species:")
    for s in sorted(missing):
        print(f"  {s}")



# Look for any EIF1-like gene names in a gastrophryne, which is missing EIF1
gff = "/project/inbreh/frog_keratin/all_ref_genomes/gffs/GCA_027917425.1_aGasCar1.pri_genomic_genomic.gff"
result = subprocess.run(
    ["grep", "-iE", "gene=eif1[^a-zA-Z0-9]|gene=EIF1[^a-zA-Z0-9]", gff],
    capture_output=True, text=True
)
print("Gastrophryne EIF1-related entries:")
print(result.stdout[:2000] if result.stdout else "  (none found)")

# Also check for any EIF1 in X. laevis with homeolog suffixes
gff2 = "/project/inbreh/frog_keratin/all_ref_genomes/gffs/GCF_017654675.1_Xenopus_laevis_v10.1_genomic_genomic.gff"
result2 = subprocess.run(
    ["awk", '-F\t', 'BEGIN{IGNORECASE=1} $3=="gene" && $9 ~ /gene=(smarce1|eif1|faim2|eif4b)/ {print $9}', gff2],
    capture_output=True, text=True
)
print("\nXenopus laevis anchor-like entries:")
print(result2.stdout[:2000] if result2.stdout else "  (none found)")





