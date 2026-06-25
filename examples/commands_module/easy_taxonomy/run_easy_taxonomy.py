"""
Easy Taxonomy — end-to-end showcase.

Classifies ~100 diverse protein sequences (bacteria, archaea, fungi,
plants, animals, viruses) against the SwissProt database, then exercises
every EasyTaxonomyParser method.

Data is downloaded automatically on first run:
  - Query: ~100 reviewed proteins from UniProt across 6 kingdoms (~50 KB)
  - Target: SwissProt via `mmseqs databases` (~250 MB download, ~1.5 GB on disk)
"""

import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt

# ── Paths (relative to this script) ─────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
DB_DIR = SCRIPT_DIR / "swissprot_db"
OUTPUT_DIR = SCRIPT_DIR / "taxonomy_output"
QUERY_FASTA = DATA_DIR / "query_diverse.fasta"
TARGET_DB = DB_DIR / "swissprot"

DATA_DIR.mkdir(exist_ok=True)
DB_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)


# ── Step 1: Download diverse query sequences from UniProt ────────
def download_queries():
    """Download ~100 reviewed proteins spanning major domains of life."""
    if QUERY_FASTA.exists() and QUERY_FASTA.stat().st_size > 0:
        n = sum(1 for line in open(QUERY_FASTA) if line.startswith(">"))
        print(f"Query FASTA already exists ({n} sequences), skipping download.")
        return

    print("Downloading query sequences from UniProt...")

    # Balanced sampling across the tree of life
    taxa = [
        ("Bacteria",      2,     40),
        ("Archaea",       2157,  15),
        ("Fungi",         4751,  15),
        ("Viridiplantae", 33090, 10),
        ("Metazoa",       33208, 10),
        ("Viruses",       10239, 10),
    ]

    all_fasta = []
    for name, taxid, size in taxa:
        query = f"taxonomy_id:{taxid} AND reviewed:true AND length:[100 TO 500]"
        params = urllib.parse.urlencode({"query": query, "format": "fasta", "size": size})
        url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
        req = urllib.request.Request(url, headers={"User-Agent": "pymmseqs-showcase/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read().decode()
                n_seqs = data.count(">")
                print(f"  {name} (taxid={taxid}): {n_seqs} sequences")
                all_fasta.append(data)
        except Exception as e:
            print(f"  Warning: failed to download {name}: {e}")

    if not all_fasta:
        print("ERROR: Could not download any query sequences.")
        sys.exit(1)

    QUERY_FASTA.write_text("".join(all_fasta))
    total = sum(chunk.count(">") for chunk in all_fasta)
    print(f"Saved {total} query sequences to {QUERY_FASTA}\n")


# ── Step 2: Download SwissProt taxonomy database ─────────────────
def setup_swissprot_db():
    """Download SwissProt with taxonomy via mmseqs databases."""
    # Check if database already exists
    db_type_file = Path(f"{TARGET_DB}.dbtype")
    taxonomy_file = Path(f"{TARGET_DB}_mapping")
    if db_type_file.exists() and taxonomy_file.exists():
        print("SwissProt database already exists, skipping download.\n")
        return

    print("Downloading SwissProt database with taxonomy (this may take a few minutes)...")
    print("  Size: ~250 MB download, ~1.5 GB on disk\n")

    tmp_dir = DB_DIR / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    subprocess.run(
        ["mmseqs", "databases", "UniProtKB/Swiss-Prot",
         str(TARGET_DB), str(tmp_dir)],
        check=True,
    )
    print("SwissProt database ready.\n")


# ── Step 3: Run easy_taxonomy ────────────────────────────────────
download_queries()
setup_swissprot_db()

from pymmseqs.commands import easy_taxonomy

print("Running easy_taxonomy...")
result = easy_taxonomy(
    fasta_file=str(QUERY_FASTA),
    target_db=str(TARGET_DB),
    tax_reports=str(OUTPUT_DIR / "tax_results"),
    s=5.7,
    lca_ranks="superkingdom,phylum,class,order,family,genus,species",
    tax_lineage=1,
)

# ── to_path ──────────────────────────────────────────────────────
print("\n=== to_path() ===")
paths = result.to_path()
for name, path in paths.items():
    print(f"  {name}: {path}")

# ── to_pandas (all 4 output types) ──────────────────────────────
print("\n=== to_pandas(output='lca') — LCA assignments ===")
lca_df = result.to_pandas(output="lca")
print(f"Shape: {lca_df.shape}")
print(lca_df.head(5).to_string())

print("\n=== to_pandas(output='report') — Kraken-style report ===")
report_df = result.to_pandas(output="report")
print(f"Shape: {report_df.shape}")
print(report_df.head(10).to_string())

print("\n=== to_pandas(output='tophit_aln') — Top-hit alignments ===")
tophit_df = result.to_pandas(output="tophit_aln")
print(f"Shape: {tophit_df.shape}")
print(tophit_df.head(5).to_string())

# ── to_list ──────────────────────────────────────────────────────
print("\n=== to_list() (first 3 entries) ===")
rows = result.to_list()
for row in rows[:3]:
    print(row)

# ── to_gen ───────────────────────────────────────────────────────
print("\n=== to_gen() (first 3 entries) ===")
gen = result.to_gen()
for _ in range(min(3, len(lca_df))):
    row = next(gen)
    print(f"  {row['query']} → taxid={row['taxid']}  {row.get('taxon_name', '')}")

# ── __len__ and __repr__ ─────────────────────────────────────────
print(f"\n=== len(result) === {len(result)}")
print(f"\n=== repr(result) ===\n{repr(result)}")

# ── summary ──────────────────────────────────────────────────────
print("\n=== summary() ===")
stats = result.summary()
for k, v in stats.items():
    print(f"  {k}: {v}")

# ── composition at multiple ranks ────────────────────────────────
for rank in ["superkingdom", "phylum", "class"]:
    print(f"\n=== composition(rank='{rank}') ===")
    comp = result.composition(rank=rank)
    print(comp.head(10).to_string())

# ── diversity at multiple ranks ──────────────────────────────────
print("\n=== diversity() ===")
for rank in ["phylum", "class", "order", "family", "genus"]:
    div = result.diversity(rank=rank)
    print(f"  {rank}: richness={div['richness']}, "
          f"shannon={div['shannon_entropy']}, "
          f"simpson={div['simpson_index']}")

# ── rank_summary ─────────────────────────────────────────────────
print("\n=== rank_summary() ===")
rank_df = result.rank_summary()
print(rank_df.to_string())

# ── unclassified_report ──────────────────────────────────────────
print("\n=== unclassified_report() ===")
uc = result.unclassified_report()
for k, v in uc.items():
    print(f"  {k}: {v}")

# ── filter_by_taxon ──────────────────────────────────────────────
print("\n=== filter_by_taxon() ===")
try:
    top_phylum = result.composition("phylum").iloc[0]["taxon"]
    filtered = result.filter_by_taxon(top_phylum, rank="phylum")
    print(f"Sequences in {top_phylum}: {len(filtered)}")
    print(filtered.head(5).to_string())
except (IndexError, ValueError) as e:
    print(f"  Skipped: {e}")

# ── to_json / to_csv ────────────────────────────────────────────
print("\n=== Exports ===")
result.to_json(OUTPUT_DIR / "lca_assignments.json")
print(f"Saved {OUTPUT_DIR / 'lca_assignments.json'}")

result.to_csv(OUTPUT_DIR / "lca_assignments.csv", output="lca")
print(f"Saved {OUTPUT_DIR / 'lca_assignments.csv'}")

result.to_csv(OUTPUT_DIR / "report.csv", output="report")
print(f"Saved {OUTPUT_DIR / 'report.csv'}")

# ── Visualizations ───────────────────────────────────────────────
print("\n=== Plots ===")

fig, ax = result.plot_composition(rank="phylum", top_n=10)
fig.savefig(OUTPUT_DIR / "composition_phylum.png", dpi=150)
print("Saved composition_phylum.png")

fig, ax = result.plot_composition(rank="class", top_n=10)
fig.savefig(OUTPUT_DIR / "composition_class.png", dpi=150)
print("Saved composition_class.png")

fig, ax = result.plot_rank_resolution()
fig.savefig(OUTPUT_DIR / "rank_resolution.png", dpi=150)
print("Saved rank_resolution.png")

fig, ax = result.plot_composition_pie(rank="phylum", top_n=8)
fig.savefig(OUTPUT_DIR / "composition_pie.png", dpi=150)
print("Saved composition_pie.png")

fig, ax = result.plot_classified_vs_unclassified()
fig.savefig(OUTPUT_DIR / "classified_vs_unclassified.png", dpi=150)
print("Saved classified_vs_unclassified.png")

fig, ax = result.plot_diversity_comparison()
fig.savefig(OUTPUT_DIR / "diversity_comparison.png", dpi=150)
print("Saved diversity_comparison.png")

plt.show()

print("\nDone — all parser methods verified successfully.")
