"""
scop_case_study.py
===================

This module implements a number of utilities to assemble training and test
datasets for benchmarking algorithms that operate on Ramachandran‐like
distributions on a sphere.  The idea is to select several sets of protein
domains from very different parts of the SCOP hierarchy (easy task),
progressively moving down to folds, superfamilies and families (moderate,
hard and challenging tasks).  For each selected domain the script can
download the corresponding PDB file and compute simple backbone bond
angles such as τ/θ (the angle between the N–Cα and C–Cα pseudo‑bonds) and,
optionally, the φ and ψ torsion angles.  Results are saved into CSV files
so that they can easily be consumed by downstream analyses.

The code does **not** require any third party Python packages – only
modules from the standard library are used.  Network access is required
to download SCOP classification files and PDB structures.  If the
SCOP files cannot be downloaded automatically (for example because of
certificate problems), they can be placed manually into the directory
specified via the ``--scop-dir`` argument.  The script will detect and
use existing files before attempting to download.

Workflow overview
-----------------

1. Download (or locate) the SCOP classification file ``dir.cla.scop.1.75.txt``.
   This file contains one record per SCOP domain and provides the
   concise classification string (sccs) and the domain identifier (sid).
   The sccs has the form ``<class>.<fold>.<superfamily>.<family>``, where
   ``class`` is a single lower‑case letter (e.g. ``a`` for all‑α proteins,
   ``b`` for all‑β proteins, etc.) and the numbers encode deeper levels
   of the hierarchy【947083788190719†L235-L245】.  By inspecting the sccs
   it is possible to group domains according to class, fold, superfamily or
   family without consulting any other files.

2. Parse the classification file into an in‑memory structure that maps
   every domain to its sccs and PDB accession.  Only classes ``a`` to
   ``g`` (the true SCOP structural classes) are considered; other
   "artificial" classes (peptides, coiled‑coil proteins, low resolution
   structures and designed proteins) are ignored【947083788190719†L246-L257】.

3. Select a number of distinct classification codes at the desired level.
   For the easy task the level is ``class`` – four different classes are
   chosen with the highest number of available domains.  For the
   moderate task the level is ``fold`` within different classes; again
   the folds with the most domains are selected.  For the hard task
   four superfamilies within the same fold are chosen.  For the
   challenging task four families within the same superfamily are selected.
   The functions ``find_top_categories`` and ``sample_domains`` implement
   this logic.  To increase reproducibility the Python random module can
   be seeded via the ``--seed`` argument.

4. For each selected domain draw up to the desired number of unique PDB
   entries (300–500 by default).  A single PDB entry may contain
   multiple domains; sampling is done at the domain level but only one
   copy of each PDB file is downloaded and used for all domains it
   contains.  The PDB coordinate files are fetched from the RCSB
   distribution at ``https://files.rcsb.org/download/{pdb_id}.pdb``.

5. Each structure is parsed using a tiny PDB parser implemented in this
   module.  For every residue the positions of the ``N``, ``CA`` and
   ``C`` atoms are extracted.  The τ/θ bond angle at ``Cα`` is then
   computed as the angle between the vectors ``N–Cα`` and ``C–Cα``.  The
   code also includes an optional torsion angle calculator for φ and
   ψ using vector algebra.

6. Aggregated results are written into CSV files organised by task and
   classification.  Each row of the CSV file contains the PDB code,
   domain id, residue number, residue name and the computed angles.

Example usage::

    python scop_case_study.py --scop-dir data/scop \
        --output-dir results --easy-size 400 --moderate-size 400 \
        --hard-size 400 --challenging-size 400

Please consult the accompanying README.md for detailed instructions.

Author: OpenAI Assistant (2025)
"""

import argparse
import csv
import os
import random
import re
import ssl
import sys
from collections import defaultdict, Counter
from dataclasses import dataclass
from math import acos, degrees, sqrt
from typing import Dict, List, Tuple, Iterable, Optional
from urllib.error import URLError
from urllib.request import urlopen


###############################################################################
# Data classes
###############################################################################

@dataclass
class DomainRecord:
    """Simple container for information about a SCOP domain.

    Attributes
    ----------
    sid : str
        The SCOP domain identifier, e.g. ``d1tpt_1``.
    sccs : str
        The concise classification string (class.fold.superfamily.family).
    pdb : str
        The 4‑letter PDB accession associated with this domain.
    chain : str
        The chain identifier within the PDB file.
    """
    sid: str
    sccs: str
    pdb: str
    chain: str


###############################################################################
# SCOP file handling
###############################################################################

SCOP_CLASSFILE_NAME = "dir.cla.scop.1.75.txt"


def download_scop_classfile(scop_dir: str, quiet: bool = False) -> str:
    """Download the SCOP classification file if it does not already exist.

    Parameters
    ----------
    scop_dir : str
        Directory in which the SCOP file should be stored.  The directory
        will be created if it does not exist.
    quiet : bool, optional
        If ``True``, suppress status messages.

    Returns
    -------
    str
        Absolute path to the downloaded (or existing) classification file.

    Notes
    -----
    The SCOP maintainers distribute parseable files for version 1.75 at
    ``https://scop.berkeley.edu/downloads/parse/``.  This function attempts
    to download ``dir.cla.scop.1.75.txt`` using urllib with SSL verification
    disabled.  If the download fails (e.g. due to certificate issues), the
    user is instructed to supply the file manually.
    """
    os.makedirs(scop_dir, exist_ok=True)
    classfile_path = os.path.join(scop_dir, SCOP_CLASSFILE_NAME)
    if os.path.exists(classfile_path):
        return classfile_path

    # Attempt to download
    url = ("https://scop.berkeley.edu/downloads/parse/" + SCOP_CLASSFILE_NAME)
    if not quiet:
        print(f"Downloading SCOP classification file from {url}...")
    # create unverified SSL context
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    try:
        with urlopen(url, context=context, timeout=30) as response:
            data = response.read()
        with open(classfile_path, "wb") as out_f:
            out_f.write(data)
        if not quiet:
            print(f"Saved SCOP classification file to {classfile_path}")
    except URLError as e:
        # Provide a helpful message and re‑raise
        msg = ("Failed to download SCOP classification file. "
               "Please download '" + SCOP_CLASSFILE_NAME + "' manually from the"
               " SCOP website (https://scop.berkeley.edu) and place it in "
               f"the directory {scop_dir}. Error: {e}")
        raise RuntimeError(msg)
    return classfile_path


def parse_scop_classfile(filename: str) -> List[DomainRecord]:
    """Parse ``dir.cla.scop.1.75.txt`` into a list of DomainRecord objects.

    The classification file contains one record per domain.  Each line has
    several whitespace‑separated fields.  The first field is the domain
    identifier (sid) prefixed by either 'd' (domain), 'e' (engineered),
    etc.  The second field is the SCOP numeric id (sunid) and the third
    field is the sccs.  The PDB accession is embedded in the sid (characters
    1–4), and the chain identifier follows the underscore.  Some sids
    correspond to engineered constructs or other non‑standard entries; those
    are ignored for the purpose of selecting natural proteins.

    Parameters
    ----------
    filename : str
        Path to the ``dir.cla.scop.1.75.txt`` file.

    Returns
    -------
    list of DomainRecord
        Parsed records.
    """
    records: List[DomainRecord] = []
    with open(filename, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # skip comments
            parts = line.split()
            # Format: sid sunid sccs residue_range
            # e.g. d1tpt_ 14919 a.1.1.2 (1-70) ...
            # Some engineered entries start with 'g' or 'd', keep all for now
            # The parseable classification files may have either three or four
            # whitespace‑separated fields before the numerical sunid/indices.  In
            # older releases the fields are: sid, sunid, sccs, residue_range.  In
            # SCOP 1.75 the third column encodes the chain (e.g. "A:") and the
            # classification string (sccs) moves to the fourth column.  Robustly
            # extract the sccs by checking for a trailing colon in parts[2].
            if len(parts) < 3:
                continue
            sid = parts[0]
            # Determine which element contains the sccs.  If the third part ends
            # with a colon (e.g. "A:" or "-"), assume the classification string
            # is the fourth element; otherwise use the third element.
            if len(parts) > 3 and parts[2].endswith(":"):
                sccs = parts[3]
            else:
                sccs = parts[2]
            # We only care about domains belonging to classes a–g (true structural classes)
            class_letter = sccs.split('.')[0].lower()
            if class_letter not in {"a", "b", "c", "d", "e", "f", "g"}:
                continue
            pdb_id = sid[1:5].lower()  # skip the initial 'd' and make lower case
            chain_id_match = re.search(r"_(\w)", sid)
            chain_id = chain_id_match.group(1) if chain_id_match else "?"
            records.append(DomainRecord(sid=sid, sccs=sccs, pdb=pdb_id, chain=chain_id))
    return records


###############################################################################
# Category selection
###############################################################################

def group_domains_by_category(records: List[DomainRecord], level: int) -> Dict[str, List[DomainRecord]]:
    """Group domains by a prefix of their sccs classification string.

    Parameters
    ----------
    records : list of DomainRecord
        List of parsed domain records.
    level : int
        Determines how many components of the sccs string are used for grouping.
        * 1 → class (e.g. 'a')
        * 2 → fold (e.g. 'a.1')
        * 3 → superfamily (e.g. 'a.1.1')
        * 4 → family (e.g. 'a.1.1.2')

    Returns
    -------
    dict
        Mapping from the category code to the list of domains falling under
        that category.
    """
    groups: Dict[str, List[DomainRecord]] = defaultdict(list)
    for rec in records:
        parts = rec.sccs.split('.')
        if len(parts) < level:
            continue
        key = '.'.join(parts[:level])
        groups[key].append(rec)
    return groups


def find_top_categories(groups: Dict[str, List[DomainRecord]], num: int, min_domains: int) -> List[str]:
    """Return the keys of the categories with the largest number of domains.

    Categories that have fewer than ``min_domains`` unique PDB entries are
    ignored.  Only the first ``num`` categories are returned.  Results are
    sorted by descending domain count.

    Parameters
    ----------
    groups : dict
        Mapping of category code to list of DomainRecord objects.
    num : int
        Number of categories to return.
    min_domains : int
        Minimum number of domains (unique PDB entries) required for a
        category to be considered.

    Returns
    -------
    list of str
        List of selected category codes.
    """
    # Compute unique PDB counts for each category
    counts: List[Tuple[str, int]] = []
    for key, recs in groups.items():
        unique_pdbs = {r.pdb for r in recs}
        count = len(unique_pdbs)
        if count >= min_domains:
            counts.append((key, count))
    # Sort by descending count
    counts.sort(key=lambda x: (-x[1], x[0]))
    return [key for key, _ in counts[:num]]


def sample_domains(records: List[DomainRecord], category_code: str, level: int, max_per_cat: int, seed: Optional[int] = None) -> List[DomainRecord]:
    """Sample up to ``max_per_cat`` unique PDB entries within a category.

    The function first filters the list of ``records`` to those whose sccs
    starts with the given ``category_code`` (e.g. 'a', 'a.1', 'a.1.1').
    It then groups those domains by PDB accession and selects at most one
    representative domain per PDB, up to ``max_per_cat``.  Sampling is
    random but reproducible if a seed is provided.

    Parameters
    ----------
    records : list of DomainRecord
        Full list of domain records.
    category_code : str
        Classification prefix specifying the category of interest.
    level : int
        Hierarchy level corresponding to ``category_code`` (1 to 4).  Only
        domains whose sccs matches exactly ``level`` components are kept.
    max_per_cat : int
        Maximum number of PDB entries to sample.
    seed : int, optional
        Seed for the random number generator.  If ``None``, sampling is
        non‑deterministic.

    Returns
    -------
    list of DomainRecord
        A list of domain records corresponding to the sampled PDB entries.
    """
    rng = random.Random(seed)
    prefix = category_code + ('' if category_code.endswith('.') else '.')
    # Filter records whose sccs starts with the prefix
    selected = [rec for rec in records if rec.sccs.startswith(prefix)]
    # Group by PDB accession
    pdb_groups: Dict[str, List[DomainRecord]] = defaultdict(list)
    for rec in selected:
        pdb_groups[rec.pdb].append(rec)
    pdb_ids = list(pdb_groups.keys())
    rng.shuffle(pdb_ids)
    chosen: List[DomainRecord] = []
    for pdb_id in pdb_ids[:max_per_cat]:
        # take the first domain for this PDB (arbitrary choice)
        chosen.append(pdb_groups[pdb_id][0])
    return chosen


###############################################################################
# PDB handling and geometry calculations
###############################################################################

def download_pdb_file(pdb_id: str, dest_dir: str) -> str:
    """Download a PDB file from the RCSB if not already present.

    Parameters
    ----------
    pdb_id : str
        Four‑letter PDB accession (lower‑case or upper‑case).
    dest_dir : str
        Directory in which to save the PDB file.

    Returns
    -------
    str
        Path to the downloaded (or existing) PDB file.
    """
    pdb_id = pdb_id.lower()
    os.makedirs(dest_dir, exist_ok=True)
    filename = os.path.join(dest_dir, f"{pdb_id}.pdb")
    if os.path.exists(filename):
        return filename
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Fetching PDB {pdb_id} from {url}")
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    try:
        with urlopen(url, context=context, timeout=30) as response:
            data = response.read()
        with open(filename, "wb") as fh:
            fh.write(data)
    except URLError as e:
        raise RuntimeError(f"Failed to download PDB {pdb_id}: {e}")
    return filename


def parse_protein_name(pdb_path: str) -> str:
    """Extract a human‐readable protein name from a PDB file.

    The function searches for `COMPND` records containing a `MOLECULE:`
    descriptor as recommended by the PDB format.  If multiple molecule
    names are present, the first is returned.  As a fallback the
    `TITLE` record is used.  If both fail the four‑letter PDB ID is
    returned.

    Parameters
    ----------
    pdb_path : str
        Path to a downloaded PDB file.

    Returns
    -------
    str
        A best‑effort description of the protein.
    """
    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0].upper()
    molecule_name: Optional[str] = None
    title_line: Optional[str] = None
    try:
        with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                record = line[0:6].strip()
                if record == 'COMPND' and 'MOLECULE:' in line:
                    # Split at 'MOLECULE:' and then at semicolon
                    parts = line.split('MOLECULE:', 1)[1].split(';', 1)
                    name = parts[0].strip()
                    if name:
                        molecule_name = name
                        break
                elif record == 'TITLE' and title_line is None:
                    title_line = line[10:].strip()
    except OSError:
        return pdb_id
    if molecule_name:
        return molecule_name
    if title_line:
        return title_line
    return pdb_id


def parse_pdb_backbone(pdb_path: str) -> List[Tuple[str, int, Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]:
    """Extract backbone atom coordinates from a PDB file.

    Only ATOM records belonging to standard residues are considered.  The
    return value is a list of tuples containing the residue name, the
    sequence number and the coordinates of the N, CA and C atoms.

    Returns
    -------
    list of tuple
        Each element is ``(res_name, res_seq, N_coords, CA_coords, C_coords)``.
    """
    backbone: Dict[Tuple[str, int], Dict[str, Tuple[float, float, float]]] = {}
    with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            # Only consider ATOM records for standard residues
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            try:
                res_seq = int(line[22:26])
            except ValueError:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            key = (chain_id, res_seq)
            if atom_name in {"N", "CA", "C"}:
                if key not in backbone:
                    backbone[key] = {"res_name": res_name}
                backbone[key][atom_name] = (x, y, z)
    # Convert to sorted list and filter residues that have all three atoms
    result: List[Tuple[str, int, Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]] = []
    for (chain_id, res_seq), atoms in sorted(backbone.items(), key=lambda x: x[0][1]):
        if all(name in atoms for name in ("N", "CA", "C")):
            result.append((atoms["res_name"], res_seq, atoms["N"], atoms["CA"], atoms["C"]))
    return result


def vector(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """Return the vector from point a to point b."""
    return (b[0] - a[0], b[1] - a[1], b[2] - a[2])


def dot(v1: Tuple[float, float, float], v2: Tuple[float, float, float]) -> float:
    """Dot product of two 3‑vectors."""
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def norm(v: Tuple[float, float, float]) -> float:
    """Euclidean norm of a 3‑vector."""
    return sqrt(dot(v, v))


def angle_between(v1: Tuple[float, float, float], v2: Tuple[float, float, float]) -> float:
    """Return the angle in degrees between two vectors."""
    n1 = norm(v1)
    n2 = norm(v2)
    if n1 == 0 or n2 == 0:
        return float('nan')
    cos_theta = max(-1.0, min(1.0, dot(v1, v2) / (n1 * n2)))
    return degrees(acos(cos_theta))


def torsion_angle(p1, p2, p3, p4) -> float:
    """Compute the torsion (dihedral) angle defined by four points.

    The torsion angle is defined by the angle between the plane formed
    by (p1, p2, p3) and the plane formed by (p2, p3, p4).  Returned
    angle is in degrees and ranges between -180 and 180.
    """
    import math
    # Vectors between the points
    b1 = vector(p1, p2)
    b2 = vector(p2, p3)
    b3 = vector(p3, p4)
    # Normal vectors to the planes
    def cross(u, v):
        return (u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0])
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    # Normalize normals
    n1n = norm(n1)
    n2n = norm(n2)
    if n1n == 0 or n2n == 0:
        return float('nan')
    n1 = (n1[0] / n1n, n1[1] / n1n, n1[2] / n1n)
    n2 = (n2[0] / n2n, n2[1] / n2n, n2[2] / n2n)
    # Compute angle between normals
    cos_phi = max(-1.0, min(1.0, dot(n1, n2)))
    phi = math.degrees(math.acos(cos_phi))
    # Determine sign using the orientation of b2 and n1 x n2
    m1 = cross(n1, n2)
    if dot(m1, b2) > 0:
        phi = -phi
    return phi


def compute_angles_for_structure(backbone: List[Tuple[str, int, Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]) -> List[Dict[str, object]]:
    """Given a list of backbone atoms, compute τ/θ, φ and ψ angles for each residue.

    Parameters
    ----------
    backbone : list
        Output of :func:`parse_pdb_backbone`.

    Returns
    -------
    list of dict
        Each dictionary contains the residue name (res_name), sequence number
        (res_seq) and the computed angles (tau, phi, psi).  The first and
        last residues in a chain will have ``nan`` torsion angles.
    """
    results = []
    n = len(backbone)
    for i, (res_name, res_seq, N_i, CA_i, C_i) in enumerate(backbone):
        # τ/θ angle at Cα: between vectors N_i→CA_i and C_i→CA_i
        v1 = vector(N_i, CA_i)
        v2 = vector(C_i, CA_i)
        tau_angle = angle_between(v1, v2)
        # Initialize torsions
        phi = float('nan')
        psi = float('nan')
        # φ: C_{i‑1}, N_i, CA_i, C_i
        if i > 0:
            prev_C = backbone[i - 1][4]
            phi = torsion_angle(prev_C, N_i, CA_i, C_i)
        # ψ: N_i, CA_i, C_i, N_{i+1}
        if i < n - 1:
            next_N = backbone[i + 1][2]
            psi = torsion_angle(N_i, CA_i, C_i, next_N)
        results.append({
            "res_name": res_name,
            "res_seq": res_seq,
            "tau": tau_angle,
            "phi": phi,
            "psi": psi
        })
    return results


###############################################################################
# Case study orchestration
###############################################################################

def run_case_study(scop_records: List[DomainRecord], output_dir: str, sizes: Dict[str, int], seed: Optional[int] = None) -> None:
    """Create datasets for each difficulty level.

    Parameters
    ----------
    scop_records : list of DomainRecord
        Parsed SCOP classification records.
    output_dir : str
        Base directory in which results and downloaded structures will be stored.
    sizes : dict
        Mapping of difficulty level ('easy', 'moderate', 'hard', 'challenging')
        to the number of PDB entries to sample per category.  Typical
        values are between 300 and 500.
    seed : int, optional
        Seed for the random number generator used when sampling.

    Notes
    -----
    The function will create one subdirectory per difficulty level.  Within
    each, a ``structures`` folder will hold the downloaded PDB files and
    a ``data.csv`` file will contain the computed angles.  An
    accompanying ``meta.txt`` file records which categories were chosen.
    """
    rng = random.Random(seed)
    os.makedirs(output_dir, exist_ok=True)
    # Group domains at different levels
    by_class = group_domains_by_category(scop_records, level=1)
    by_fold = group_domains_by_category(scop_records, level=2)
    by_superfamily = group_domains_by_category(scop_records, level=3)
    by_family = group_domains_by_category(scop_records, level=4)
    # Easy: pick four classes
    easy_classes = find_top_categories(by_class, num=4, min_domains=sizes['easy'])
    # Moderate: pick four folds from different classes
    # For fairness we avoid selecting folds whose class overlaps with others
    moderate_folds = []
    used_classes: set = set()
    # Sort all folds by domain count
    sorted_folds = sorted([(k, len({r.pdb for r in v})) for k, v in by_fold.items()], key=lambda x: -x[1])
    for code, count in sorted_folds:
        class_letter = code.split('.')[0]
        if class_letter in used_classes:
            continue
        if count >= sizes['moderate']:
            moderate_folds.append(code)
            used_classes.add(class_letter)
        if len(moderate_folds) == 4:
            break
    # Hard: pick four superfamilies from within the first selected fold
    hard_superfamilies = []
    if moderate_folds:
        fold_for_hard = moderate_folds[0]
        # Collect superfamilies belonging to this fold
        sf_groups = {k: v for k, v in by_superfamily.items() if k.startswith(fold_for_hard + '.')}
        hard_superfamilies = find_top_categories(sf_groups, num=4, min_domains=sizes['hard'])
    # Challenging: pick four families from within the first selected superfamily
    challenging_families = []
    if hard_superfamilies:
        sf_for_chal = hard_superfamilies[0]
        fam_groups = {k: v for k, v in by_family.items() if k.startswith(sf_for_chal + '.')}
        challenging_families = find_top_categories(fam_groups, num=4, min_domains=sizes['challenging'])
    # Summary of chosen categories
    selection = {
        'easy': easy_classes,
        'moderate': moderate_folds,
        'hard': hard_superfamilies,
        'challenging': challenging_families
    }
    # Process each difficulty level
    for level_name, categories in selection.items():
        if not categories:
            print(f"Warning: no categories selected for {level_name}")
            continue
        level_dir = os.path.join(output_dir, level_name)
        struct_dir = os.path.join(level_dir, 'structures')
        os.makedirs(struct_dir, exist_ok=True)
        # Write meta information
        with open(os.path.join(level_dir, 'meta.txt'), 'w') as meta_f:
            meta_f.write(f"Selected categories for {level_name}:\n")
            for cat in categories:
                meta_f.write(f"  {cat}\n")
        # Prepare CSV output
        csv_path = os.path.join(level_dir, 'data.csv')
        with open(csv_path, 'w', newline='') as csvfile:
            # Additional fields: class_letter and protein_name
            fieldnames = [
                'pdb_id', 'protein_name', 'domain_id', 'class', 'category',
                'res_name', 'res_seq', 'tau', 'phi', 'psi'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for cat in categories:
                # Determine hierarchy level based on number of components
                level = cat.count('.') + 1
                records_for_cat = sample_domains(
                    scop_records, cat, level, sizes[level_name], seed=rng.randint(0, 2**32 - 1)
                )
                for rec in records_for_cat:
                    pdb_path = download_pdb_file(rec.pdb, struct_dir)
                    protein_name = parse_protein_name(pdb_path)
                    backbone = parse_pdb_backbone(pdb_path)
                    angles = compute_angles_for_structure(backbone)
                    class_letter = rec.sccs.split('.')[0]
                    for angle in angles:
                        writer.writerow({
                            'pdb_id': rec.pdb,
                            'protein_name': protein_name,
                            'domain_id': rec.sid,
                            'class': class_letter,
                            'category': cat,
                            'res_name': angle['res_name'],
                            'res_seq': angle['res_seq'],
                            'tau': angle['tau'],
                            'phi': angle['phi'],
                            'psi': angle['psi']
                        })


###############################################################################
# Command line interface
###############################################################################

def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Generate SCOP case study datasets")
    parser.add_argument('--scop-dir', default='scop_data', help='Directory for SCOP parse files')
    parser.add_argument('--output-dir', default='results', help='Directory to write the datasets and downloaded structures')
    parser.add_argument('--easy-size', type=int, default=400, help='Number of PDB entries to sample per class for the easy task')
    parser.add_argument('--moderate-size', type=int, default=400, help='Number of PDB entries to sample per fold for the moderate task')
    parser.add_argument('--hard-size', type=int, default=400, help='Number of PDB entries to sample per superfamily for the hard task')
    parser.add_argument('--challenging-size', type=int, default=400, help='Number of PDB entries to sample per family for the challenging task')
    parser.add_argument('--seed', type=int, default=None, help='Random seed for reproducibility')
    args = parser.parse_args(argv)
    # Download or locate the SCOP classification file
    try:
        classfile_path = download_scop_classfile(args.scop_dir)
    except RuntimeError as e:
        print(e, file=sys.stderr)
        return 1
    # Parse the classification file
    print(f"Parsing SCOP classification from {classfile_path}...")
    records = parse_scop_classfile(classfile_path)
    print(f"Parsed {len(records)} domain records")
    # Prepare sizes dict
    sizes = {
        'easy': args.easy_size,
        'moderate': args.moderate_size,
        'hard': args.hard_size,
        'challenging': args.challenging_size
    }
    run_case_study(records, args.output_dir, sizes, seed=args.seed)
    print("Case study generation complete.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
