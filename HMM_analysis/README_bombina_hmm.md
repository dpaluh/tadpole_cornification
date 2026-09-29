HMM search for keratins in Bombina


Get human and mouse alpha-keratins from Type-I and Type-II families from UniProt Swiss-Prot. Query by protein family field.


Get out metadata tsv:

```
WORKDIR=/project/inbreh/frog_keratin/seed_sequences
mkdir -p $WORKDIR
cd $WORKDIR

BASE_URL="https://rest.uniprot.org/uniprotkb/stream"

# Fields to retrieve in TSV for inspection
FIELDS="accession,gene_names,protein_name,organism_name,length,cc_subcellular_location"

for SPECIES in "9606" "10090"; do
    if [ "$SPECIES" = "9606" ]; then
        LABEL="human"
    else
        LABEL="mouse"
    fi

    echo "Fetching $LABEL metadata..."

    # Type I
    curl -s "${BASE_URL}?query=(organism_id:${SPECIES}+AND+reviewed:true+AND+protein_name:%22Keratin,+type+I%22)&fields=${FIELDS}&format=tsv" \
        > ${LABEL}_typeI_metadata.tsv

    # Type II
    curl -s "${BASE_URL}?query=(organism_id:${SPECIES}+AND+reviewed:true+AND+protein_name:%22Keratin,+type+II%22)&fields=${FIELDS}&format=tsv" \
        > ${LABEL}_typeII_metadata.tsv

    echo "$LABEL type I entries:  $(tail -n +2 ${LABEL}_typeI_metadata.tsv | wc -l)"
    echo "$LABEL type II entries: $(tail -n +2 ${LABEL}_typeII_metadata.tsv | wc -l)"
done
```


Download:

```
BASE_URL="https://rest.uniprot.org/uniprotkb/stream"

for SPECIES in "9606" "10090"; do
    if [ "$SPECIES" = "9606" ]; then LABEL="human"; else LABEL="mouse"; fi

    curl -s "${BASE_URL}?query=(organism_id:${SPECIES}+AND+reviewed:true+AND+protein_name:%22Keratin,+type+I%22)&format=fasta" \
        > ${LABEL}_typeI.fasta

    curl -s "${BASE_URL}?query=(organism_id:${SPECIES}+AND+reviewed:true+AND+protein_name:%22Keratin,+type+II%22)&format=fasta" \
        > ${LABEL}_typeII.fasta

    echo "$LABEL type I:  $(grep -c '^>' ${LABEL}_typeI.fasta)"
    echo "$LABEL type II: $(grep -c '^>' ${LABEL}_typeII.fasta)"
done
```

Remove krtap from mouse:


```
# Confirm the contaminant is there
grep "O08884" mouse_typeII.fasta

# Remove it
awk '/^>/{keep=1} /O08884/{keep=0} keep{print}' mouse_typeII.fasta > mouse_typeII_clean.fasta

# Verify count dropped by 1
grep -c '^>' mouse_typeII_clean.fasta
```


Concatenate all seed proteins for each:

```
cat human_typeI.fasta mouse_typeI.fasta          > seeds_typeI.fasta
cat human_typeII.fasta mouse_typeII_clean.fasta  > seeds_typeII.fasta

echo "Type I seeds:  $(grep -c '^>' seeds_typeI.fasta)"
echo "Type II seeds: $(grep -c '^>' seeds_typeII.fasta)"
```



Align each cluster with MAFFT:


```
# run this from an interactive session out of good practice, should be lightweight
cd /project/inbreh/frog_keratin/seed_sequences

salloc -A inbreh -t 0-03:00 --mem=16G --cpus-per-task=4

#load up mafft module
module load mafft/7.526

mafft --auto --thread 4 seeds_typeI.fasta  > seeds_typeI_aln.fasta
mafft --auto --thread 4 seeds_typeII.fasta > seeds_typeII_aln.fasta


# load conda to use trimal to clean out columns with too many gaps:

module load miniconda3/24.3.0

# create environment and install trimal, only once, uncomment if needed again
# conda create -n trimal bioconda::trimal -y
conda activate trimal

# Run trimal on each alignment:
trimal -in seeds_typeI_aln.fasta  -out seeds_typeI_trim.fasta  -gt 0.2
trimal -in seeds_typeII_aln.fasta -out seeds_typeII_trim.fasta -gt 0.2

# Build hmm for each krt family:
hmmbuild seeds_typeI.hmm  seeds_typeI_trim.fasta
hmmbuild seeds_typeII.hmm seeds_typeII_trim.fasta

# press each hmm:
hmmpress seeds_typeI.hmm
hmmpress seeds_typeII.hmm
```



Run hmm search against Bombina genome and human and mouse as a check:

```
PROTDIR=/project/inbreh/frog_keratin/all_ref_genomes/protein_fas
OUTDIR=/project/inbreh/frog_keratin/hmm_searches
HMMDIR=/project/inbreh/frog_keratin/seed_sequences

mkdir -p $OUTDIR

for FAA in $PROTDIR/*.faa; do
    BASENAME=$(basename $FAA _protein.faa)
    echo "Searching $BASENAME..."
    hmmsearch --cpu 4 -E 1e-5 \
        --tblout ${OUTDIR}/${BASENAME}_typeI.tbl \
        $HMMDIR/seeds_typeI.hmm $FAA

    hmmsearch --cpu 4 -E 1e-5 \
        --tblout ${OUTDIR}/${BASENAME}_typeII.tbl \
        $HMMDIR/seeds_typeII.hmm $FAA
done

```




- `parse_hmm_out.py` parses the output from the hmm searches - run it interactively from vscode: 
	- Results from this are in `/project/inbreh/frog_keratin/hmm_processed/`
	- Number of mouse type II keratins is less than reported other places, because krt74 is annotated as a transcribed pseudogene, not a functional protein, so not in protein fasta
	- We are seeing krt18 consistently as a typeI in the type II chromosome, as expected
	- others that are scattered onto different chromosomes seem to be pairs of krt18 & krt 8, which form a dimer


- `fasta_from_hmmout.py` gets those sequences and writes them to a fasta file per keratin type. 

