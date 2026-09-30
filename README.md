# NAU Comparative Genomics INF 515

## Command Line BLAST Lab

**Prof. Marc Tollis**  
marc.tollis@nau.edu

## Goal

The goal of this lab is to use command-line BLAST, rather than web BLAST, to compare mouse and zebrafish protein sequences and explore how sequence similarity can be used to identify candidate homologs and putative orthologs.

By the end of the lab, you should be able to:

- build a local BLAST database;
- run `blastp` from the command line;
- interpret standard and tabular BLAST output;
- compare percent identity, alignment length, E value, and bit score;
- identify the best zebrafish hit for each mouse protein;
- perform a reciprocal best-hit search as a simple heuristic for putative orthology.

---

## 1. Set up the project directory

```bash
export PROJ_DIR=/scratch/mt2245/CompGenomicsCourse/fall26/blast
mkdir $PROJ_DIR
cd $PROJ_DIR/
```

Download the two protein FASTA files from this GitHub repository:

```bash
wget https://github.com/marctollis/INF-515-BLAST-Lab/raw/main/mouse.1.protein.faa.gz
wget https://github.com/marctollis/INF-515-BLAST-Lab/raw/main/zebrafish.1.protein.faa.gz
```

The repository should contain only:

```text
README.md
mouse.1.protein.faa.gz
zebrafish.1.protein.faa.gz
```

Look at the files in the directory:

```bash
ls -l
```

Unzip the data:

```bash
gunzip *.faa.gz
```

Inspect the beginning of the mouse protein file:

```bash
head mouse.1.protein.faa
```

### Questions

- What file format is this?
- What kind of sequence data does it contain: nucleotide or protein?
- How can you tell what the sequence type is from the filename?

---

## 2. Start with a small query set

Take the first few mouse sequences and save them to a new file:

```bash
head -n 11 mouse.1.protein.faa > mm-first.faa
```

---

## 3. Let's BLAST

Load the BLAST module:

```bash
module load blast
```

Format the zebrafish protein FASTA as a BLAST database:

```bash
srun makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

Run `blastp` using the mouse sequences as queries:

```bash
srun blastp -query mm-first.faa -db zebrafish.1.protein.faa
```

This produces a large amount of output on the screen. Redirect it to a file instead:

```bash
srun blastp \
-query mm-first.faa \
-db zebrafish.1.protein.faa \
-out mm-first.x.zebrafish.txt
```

View the output:

```bash
less mm-first.x.zebrafish.txt
```

### Question

Why is standard BLAST output difficult to parse when we want to compare many sequences?

---

## 4. Run BLAST with a larger query set

Take a larger set of mouse proteins:

```bash
head -498 mouse.1.protein.faa > mm-second.fa
```

Run BLAST and save the standard output:

```bash
srun blastp \
-query mm-second.fa \
-db zebrafish.1.protein.faa \
-out mm-second.x.zebrafish.txt
```

Now run the same search using BLAST tabular output (`-outfmt 6`):

```bash
srun blastp \
-query mm-second.fa \
-db zebrafish.1.protein.faa \
-out mm-second.x.zebrafish.tsv \
-outfmt 6
```

The default `outfmt 6` columns are:

1. query sequence ID
2. subject sequence ID
3. percent identity
4. alignment length
5. mismatches
6. gap openings
7. query start
8. query end
9. subject start
10. subject end
11. E value
12. bit score

---

## 5. Post-BLAST analysis

How many matches are there in total?

```bash
wc -l mm-second.x.zebrafish.tsv
```

Sort the results by bit score (column 12) and inspect the highest- and lowest-scoring hits:

```bash
sort -rn -k12 mm-second.x.zebrafish.tsv | head -25
sort -rn -k12 mm-second.x.zebrafish.tsv | tail -25
```

The `r` option sorts in reverse order and `n` forces a numeric sort.

### Questions

- What is the relationship between bit score (column 12) and E value (column 11)?
- What does a high bit score mean?
- What does a low E value mean?

---

## 6. Filter by percent identity

Collect hits with at least 65% protein identity and sort them by bit score:

```bash
awk '$3 >= 65.00' mm-second.x.zebrafish.tsv | sort -rn -k12
```

Count them:

```bash
awk '$3 >= 65.00' mm-second.x.zebrafish.tsv | wc -l
```

### Questions

- What is the relationship between E value, bit score, and percent identity?
- Why can a hit with 65% or greater protein identity still have a relatively low bit score and high E value?
- Check column 4: what role does alignment length play?

**Important:** percent identity alone is not enough to evaluate the strength of a BLAST match.

---

# Reciprocal Best Hit Extension

So far we have identified zebrafish proteins that are similar to mouse proteins.

A high-scoring BLAST hit does **not** automatically prove that two genes are orthologs. Paralogs can also be highly similar.

One simple heuristic for identifying putative orthologs is the **reciprocal best hit (RBH)**:

**mouse protein → best zebrafish hit → BLAST back against mouse**

If the zebrafish protein's best hit back in mouse is the original mouse protein, the pair is a reciprocal best hit.

---

## 7. Find the best zebrafish hit for each mouse query

Sort first by mouse query ID and then by decreasing bit score:

```bash
sort -k1,1 -k12,12nr mm-second.x.zebrafish.tsv > sorted.tsv
```

Keep only the first, highest-bit-score hit for each mouse query:

```bash
awk '!seen[$1]++' sorted.tsv > best_hits.tsv
```

Count them:

```bash
wc -l best_hits.tsv
```

Inspect the first few:

```bash
head best_hits.tsv
```

In `best_hits.tsv`:

- column 1 = mouse query;
- column 2 = best zebrafish subject;
- column 11 = E value;
- column 12 = bit score.

---

## 8. How many unique zebrafish proteins are represented?

Multiple mouse proteins may have the same zebrafish protein as their best hit.

```bash
awk '{print $2}' best_hits.tsv | sort -u | wc -l
```

### Discussion

Why might many mouse proteins point to the same zebrafish best hit?

Possible biological explanations include:

- gene duplication;
- paralogy;
- one-to-many or many-to-one relationships;
- lineage-specific gene loss;
- different rates of sequence evolution.

---

## 9. Extract the zebrafish best-hit sequences

Save the unique zebrafish IDs:

```bash
awk '{print $2}' best_hits.tsv | sort -u > zebrafish-besthit-ids.txt
```

Check them:

```bash
head zebrafish-besthit-ids.txt
```

Rebuild the zebrafish BLAST database with sequence IDs indexed:

```bash
srun makeblastdb \
-in zebrafish.1.protein.faa \
-dbtype prot \
-parse_seqids
```

Extract those proteins from the zebrafish database:

```bash
blastdbcmd \
-db zebrafish.1.protein.faa \
-entry_batch zebrafish-besthit-ids.txt \
> zebrafish-besthits.faa
```

Inspect the resulting FASTA:

```bash
head zebrafish-besthits.faa
```

---

## 10. BLAST the zebrafish proteins back against mouse

Build the mouse BLAST database with indexed sequence IDs:

```bash
srun makeblastdb \
-in mouse.1.protein.faa \
-dbtype prot \
-parse_seqids
```

Run the reverse BLAST search:

```bash
srun blastp \
-query zebrafish-besthits.faa \
-db mouse.1.protein.faa \
-out zebrafish.x.mouse.tsv \
-outfmt 6
```

Now we have:

- **forward search:** mouse → zebrafish
- **reverse search:** zebrafish → mouse

---

## 11. Find the best mouse hit for each zebrafish protein

```bash
sort -k1,1 -k12,12nr zebrafish.x.mouse.tsv \
| awk '!seen[$1]++' \
> reverse_best_hits.tsv
```

Count the reverse best hits:

```bash
wc -l reverse_best_hits.tsv
```

`reverse_best_hits.tsv` contains:

**zebrafish query → best mouse hit**

---

## 12. Identify reciprocal best-hit pairs

Make a simple two-column table of the forward best-hit pairs:

```bash
cut -f1,2 best_hits.tsv | sort -u > forward_pairs.tsv
```

This has the form:

```text
mouse_protein    zebrafish_protein
```

Now make a two-column table from the reverse search and flip the columns so that it has the same orientation:

```bash
awk '{print $2 "\t" $1}' reverse_best_hits.tsv \
| sort -u \
> reverse_pairs_flipped.tsv
```

This also has the form:

```text
mouse_protein    zebrafish_protein
```

Find pairs present in both files:

```bash
comm -12 forward_pairs.tsv reverse_pairs_flipped.tsv \
> reciprocal_best_hits.tsv
```

Count them:

```bash
wc -l reciprocal_best_hits.tsv
```

View them:

```bash
cat reciprocal_best_hits.tsv
```

---

## 13. Discussion: what did we learn?

In our example, many mouse proteins had strong zebrafish BLAST hits, but very few pairs satisfied the strict reciprocal-best-hit criterion.

### Questions

- Why can two proteins be highly similar without being orthologs?
- Why can multiple mouse proteins have the same zebrafish best hit?
- Why might a true orthologous relationship fail the strict RBH test?
- How can gene duplication and gene loss complicate pairwise BLAST-based orthology inference?
- Why might genome-scale tools use orthogroups and gene trees rather than only reciprocal BLAST hits?

---

## Main takeaway

**BLAST → sequence similarity / candidate homologs**

**Reciprocal best hit → simple heuristic for putative orthology**

**Orthogroups + gene trees → more complete inference of orthology and paralogy**

A useful next question is:

> Why might we use a tool such as OrthoFinder rather than relying only on reciprocal BLAST hits?
