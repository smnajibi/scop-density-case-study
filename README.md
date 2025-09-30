# SCOP Density Case Study

This repository contains Python code to assemble benchmarking datasets
for evaluating algorithms that estimate Ramachandran‐plot densities on
the surface of a sphere.  The central idea is to select protein domains
from different parts of the SCOP hierarchical classification and
compare how an algorithm performs when the training data are
increasingly fine grained.  The tasks range from **easy** (very distinct
structural classes) to **challenging** (families within a single
superfamily).

**Important:** the code does **not** bundle any SCOP or PDB data.  It
expects the user to either place the SCOP parseable file
`dir.cla.scop.1.75.txt` into the `scop_data/` directory or allow the
script to download it automatically.  Likewise, PDB files are fetched
on demand from the RCSB archive.  If network access is restricted,
download the SCOP file and required PDB files manually and place them
into the appropriate directories before running the script.

## Background

The Structural Classification of Proteins (SCOP) groups protein
domains into a hierarchy based on structural and evolutionary
relationships.  At the top of the tree are **classes** (all‐α,
all‐β, α/β etc.), followed by **folds**, **superfamilies** and
**families**.  Classes roughly reflect secondary structure content,
while deeper levels capture increasingly close evolutionary
relationships【947083788190719†L235-L303】.  The SCOP release 1.75 is the
latest fully curated version; its parseable files are distributed by
the SCOPe project.  This repository uses the `dir.cla.scop.1.75.txt`
file to infer the classification of each domain.

## Overview of the workflow

1. **Classification parsing:** The script reads the SCOP
   classification file (`dir.cla.scop.1.75.txt`) and extracts for each
   domain the **SCOP concise classification string** (sccs) and the
   associated 4‐letter PDB code.  Only domains belonging to true
   structural classes (`a` through `g`) are retained【947083788190719†L246-L257】.

2. **Category selection:** For each difficulty level the code
   identifies suitable SCOP categories:

   - **Easy:** four distinct classes with the largest numbers of
     available domains (for example, all‐α, all‐β, α/β (parallel) and
     α+β (antiparallel))【947083788190719†L235-L245】.
   - **Moderate:** four folds from different classes.  The folds with
     the most unique PDB entries are chosen.
   - **Hard:** four superfamilies within the first chosen fold.
   - **Challenging:** four families within the first chosen
     superfamily.

3. **Sampling:** For each selected category the script samples 300–500
   unique PDB entries (by default) and uses only one domain per PDB
   entry.  The random seed can be set to ensure reproducibility.

4. **Downloading PDB files:** PDB coordinate files are downloaded from
   `https://files.rcsb.org/download/{pdb_id}.pdb` and stored in
   subdirectories under the output directory.  If a file is already
   present it is reused.

5. **Computing backbone angles:** The parser extracts the N, Cα and C
   atoms for each residue and computes the **τ/θ** angle (between
   vectors N→Cα and C→Cα).  It also calculates the **φ** and **ψ**
   torsion angles via simple vector algebra.  The implementation is
   self–contained and does not rely on external packages.

6. **Exporting results:** For each difficulty level a `data.csv`
   containing the angles for all sampled structures is written, along
   with a `meta.txt` listing the selected SCOP categories.

## Installation and prerequisites

The script is intentionally lightweight and depends only on the Python
standard library.  It runs under Python 3.7 or later.  To obtain the
SCOP classification file you need network access, or you can download
`dir.cla.scop.1.75.txt` from the SCOPe server and place it into the
`scop_data/` directory.  The PDB files are downloaded from the RCSB
archive.

## Running the case study generation

To generate datasets with 400 PDB entries per category (adjust the
sizes as desired):

```bash
python scop_case_study.py \
    --scop-dir scop_data \
    --output-dir results \
    --easy-size 400 \
    --moderate-size 400 \
    --hard-size 400 \
    --challenging-size 400 \
    --seed 12345
```

This command will:

1. Download `dir.cla.scop.1.75.txt` into `scop_data/` if it is not already there.
2. Parse the classification and identify suitable categories for the four
   difficulty levels.
3. Create subdirectories under `results/` named `easy/`,
   `moderate/`, `hard/` and `challenging/`.
4. Download the required PDB files into `structures/` subfolders.
5. Compute bond and torsion angles and write them into CSV files.

If you encounter download issues (e.g. due to certificate errors),
obtain the SCOP file and PDB files manually and place them into the
specified directories.  The script will reuse existing files.

## License

The code in this repository is released under the MIT License.  The
SCOP classification files and PDB data are subject to their own
licenses – please consult the respective providers for details.