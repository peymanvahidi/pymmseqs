"""
Easy RBH (Reciprocal Best Hits) — end-to-end showcase.

Finds orthologs between 500 human and 500 rat SwissProt proteins,
then exercises every EasyRbhParser method to verify the pipeline.
"""

# ── Run easy_rbh ─────────────────────────────────────────────────
from pymmseqs.commands import easy_rbh

result = easy_rbh(
    query_fasta="rat.fasta",
    target_fasta_or_db="human.fasta",
    alignment_file="rbh_output/rbh_results",
    s=7.5,                # sensitive search
    min_seq_id=0.3,       # keep orthologs with ≥30% identity
    format_mode=4,        # BLAST-TAB with column headers
)

# ── to_path ──────────────────────────────────────────────────────
print("\n=== to_path() ===")
print(result.to_path())

# ── to_pandas ────────────────────────────────────────────────────
print("\n=== to_pandas() ===")
df = result.to_pandas()
print(f"Shape: {df.shape}")
print(df.head(10).to_string())

# ── to_list ──────────────────────────────────────────────────────
print("\n=== to_list() (first 3 entries) ===")
rows = result.to_list()
for row in rows[:3]:
    print(row)

# ── to_gen ───────────────────────────────────────────────────────
print("\n=== to_gen() (first 3 entries) ===")
gen = result.to_gen()
for _ in range(min(3, len(df))):
    row = next(gen)
    print(f"  {row['query']} → {row['target']}  identity={row['fident']}")

# ── summary ──────────────────────────────────────────────────────
print("\n=== summary() ===")
stats = result.summary()
for k, v in stats.items():
    print(f"  {k}: {v}")

# ── to_json ──────────────────────────────────────────────────────
print("\n=== to_json() ===")
result.to_json("rbh_output/rbh_results.json")
print("Saved rbh_output/rbh_results.json")

# ── to_csv ───────────────────────────────────────────────────────
print("\n=== to_csv() ===")
result.to_csv("rbh_output/rbh_results.csv")
print("Saved rbh_output/rbh_results.csv")

# ── Visualizations ───────────────────────────────────────────────
import matplotlib.pyplot as plt
from pathlib import Path

output_dir = Path(result.to_path()).parent

print("\n=== Plots ===")

fig, ax = result.plot_identity_distribution()
fig.savefig(output_dir / "identity_distribution.png", dpi=150)
print("Saved identity_distribution.png")

fig, ax = result.plot_evalue_distribution()
fig.savefig(output_dir / "evalue_distribution.png", dpi=150)
print("Saved evalue_distribution.png")

fig, ax = result.plot_alignment_length()
fig.savefig(output_dir / "alignment_length.png", dpi=150)
print("Saved alignment_length.png")

fig, ax = result.plot_score_distribution(score="bits")
fig.savefig(output_dir / "bits_distribution.png", dpi=150)
print("Saved bits_distribution.png")

plt.show()

print("\nDone — all parser methods verified successfully.")
