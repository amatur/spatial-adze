import tskit
import pyslim
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Load SLiM .trees ===
ts = tskit.load("slim/out/neutral_100000.trees")
print("Loaded tree sequence with", ts.num_individuals, "individuals")

# === Overlay mutations ===
mutation_rate = 1e-7
next_id = pyslim.next_slim_mutation_id(ts)
mutated_ts = msprime.sim_mutations(
    ts, rate=mutation_rate,
    model=msprime.SLiMMutationModel(type=0, next_id=next_id),
    keep=True
)
print("Mutations overlaid with rate =", mutation_rate)

# === Sample individuals ===
sample_size = 2000
all_inds = [ind.id for ind in mutated_ts.individuals() if len(ind.nodes) > 0]
sampled_inds = np.random.choice(all_inds, min(sample_size, len(all_inds)), replace=False)
print(f"Sampled {len(sampled_inds)} individuals.")

# === Select 10 variants closest to genome center ===
selected_variants = []
genome_middle = mutated_ts.sequence_length / 2

for v in mutated_ts.variants():
    dist = abs(v.site.position - genome_middle)
    selected_variants.append((dist, v))
    if len(selected_variants) >= 1000:  # just collect a subset to sort from
        break

# Get 10 closest to the middle
NUMBER_OF_VARIANTS = 5
selected_variants.sort(key=lambda x: x[0])
selected_variants = [v for _, v in selected_variants[:NUMBER_OF_VARIANTS]]

print(f"Selected {len(selected_variants)} central variants for haplotypes.")

# === Build haplotypes safely (no genotype_matrix) ===
haplo_strings = []

for ind_id in sampled_inds:
    ind = mutated_ts.individual(ind_id)
    haplo = ""
    for node in ind.nodes:
        hap = "".join(str(v.genotypes[node]) for v in selected_variants)
        haplo += hap  # Concatenate both haploid genomes
    haplo_strings.append(haplo)

print("Extracted haplotypes.")

# === Assign haplotype clusters ===
haplo_to_cluster = {h: i for i, h in enumerate(sorted(set(haplo_strings)))}
haplo_clusters = [haplo_to_cluster[h] for h in haplo_strings]

print(f"Unique haplotypes found: {len(haplo_to_cluster)}")

# === Build DataFrame for plotting ===
data = []
for ind_id, cluster in zip(sampled_inds, haplo_clusters):
    ind = mutated_ts.individual(ind_id)
    x, y = ind.location[:2]
    data.append([x, y, cluster])

df = pd.DataFrame(data, columns=["x", "y", "haplotype_cluster"])

# === Plot ===
plt.figure(figsize=(8, 8))
scatter = plt.scatter(df["x"], df["y"], c=df["haplotype_cluster"], cmap="tab20", s=20)
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.title("Individuals Colored by Haplotype Cluster (10 Central Variants)")
plt.colorbar(scatter, label="Haplotype Cluster")
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig("haplotype_cluster_plot.png", dpi=300)
print("Plot saved as haplotype_cluster_plot.png")
plt.show()


# Save site positions used
site_positions = [v.site.position for v in selected_variants]

# Add haplotype strings to DataFrame
df["haplotype"] = haplo_strings

# Save to CSV
csv_filename = "haplotype_clusters.csv"
df.to_csv(csv_filename, index=False)

print(f"Saved data to {csv_filename}")
print(f"Site positions used for haplotypes: {site_positions}")
