import tskit
import pyslim
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Step 1: Load Tree Sequence with Mutations ===
input_ts = "slim/out/neutral_1500.trees"  # Replace with your .trees file
ts = tskit.load(input_ts)

print("Loaded tree sequence with", ts.num_individuals, "individuals")
# Overlay Mutations
mutation_rate = 1e-7
next_id = pyslim.next_slim_mutation_id(ts)
ts = msprime.sim_mutations(ts, rate=mutation_rate,
                                   model=msprime.SLiMMutationModel(type=0, next_id=next_id),
                                   keep=True)

print("Mutations overlaid with rate =", mutation_rate)


# === Step 2: Uniformly Sample 1000 Individuals ===
sample_size = 100
all_inds = [ind.id for ind in ts.individuals() if len(ind.nodes) > 0]
n_sample = min(sample_size, len(all_inds))
sampled_inds = np.random.choice(all_inds, n_sample, replace=False)

print("Sampled", n_sample, "individuals.")

# === Step 3: Extract Haplotypes for Each Individual ===
# Extract the genotype matrix (rows = sites, columns = nodes)
G = ts.genotype_matrix()

haplo_strings = []

# Get first 2 segregating sites
sites = []
for v in ts.variants():
    sites.append(v.genotypes)
    if len(sites) == 2:
        break

if len(sites) < 2:
    print("WARNING: Less than 2 segregating sites!")
    
for ind_id in sampled_inds:
    ind = ts.individual(ind_id)
    haplo = ""
    for node in ind.nodes:
        geno = "".join([str(site[node]) for site in sites])
        haplo += geno  # Concatenate for diploid
    haplo_strings.append(haplo)

# for ind_id in sampled_inds:
#     ind = ts.individual(ind_id)
#     haplo = ""
#     for node in ind.nodes:  # For both genomes of diploid
#         node_genotype = G[:, node]  # All site genotypes for this node
#         #haplo += "".join(node_genotype.astype(str))  # Convert to string
#         haplo += "".join(node_genotype[:2].astype(str))  # Only first 2 sites

#     haplo_strings.append(haplo)


print("Unique haplotypes found:", len(set(haplo_strings)))

# === Step 4: Map Haplotypes to Cluster IDs ===
haplo_to_cluster = {h: i for i, h in enumerate(set(haplo_strings))}
haplo_clusters = [haplo_to_cluster[h] for h in haplo_strings]

# === Step 5: Extract Positions and Build DataFrame ===
data = []
for ind_id, cluster in zip(sampled_inds, haplo_clusters):
    ind = ts.individual(ind_id)
    pos = ind.location[:2]  # x, y coordinates
    data.append([pos[0], pos[1], cluster])

df = pd.DataFrame(data, columns=["x", "y", "haplotype_cluster"])

# === Step 6: Plot ===
plt.figure(figsize=(8, 8))
scatter = plt.scatter(
    df["x"], df["y"],
    c=df["haplotype_cluster"], cmap="tab20", s=20
)

plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.title("Individuals Colored by Haplotype Cluster")
plt.colorbar(scatter, label="Haplotype Cluster")
plt.gca().set_aspect('equal')
plt.tight_layout()

# Save the plot to a PNG file
output_plot = "haplotype_clusters.png"
plt.savefig(output_plot, dpi=300)
print("Plot saved as", output_plot)

plt.show()
