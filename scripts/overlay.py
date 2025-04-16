import tskit
import pyslim
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load .trees from SLiM
ts = tskit.load("slim/out/neutral_100000.trees")
print("Loaded tree sequence with", ts.num_individuals, "individuals")

# Overlay Mutations
mutation_rate = 1e-7
next_id = pyslim.next_slim_mutation_id(ts)
mutated_ts = msprime.sim_mutations(ts, rate=mutation_rate,
                                   model=msprime.SLiMMutationModel(type=0, next_id=next_id),
                                   keep=True)

# mutated_ts = msprime.sim_mutations(ts, rate=mutation_rate,
#                                    model=msprime.InfiniteAlleles(),
#                                    keep=True)


print("Mutations overlaid with rate =", mutation_rate)

# mutated_ts.dump("output_100000_mutated.trees")

# Sample 1000 individuals uniformly
sample_size = 1000
all_inds = [ind.id for ind in mutated_ts.individuals() if len(ind.nodes) > 0]
n_sample = min(sample_size, len(all_inds))
sampled_inds = np.random.choice(all_inds, n_sample, replace=False)

print("Sampled", n_sample, "individuals.")

# Extract positions + genotypes
variant = next(mutated_ts.variants())  # First segregating site


# # Get all variant positions
# positions = [v.site.position for v in mutated_ts.variants()]

# # Find site closest to genome center
# genome_middle = mutated_ts.sequence_length / 2
# middle_site_idx = np.argmin([abs(p - genome_middle) for p in positions])

# print(f"Selected site position: {positions[middle_site_idx]}")

# # Extract the variant at this site
# variant = list(mutated_ts.variants())[middle_site_idx]

# Pick the variant closest to the middle of the genome
genome_middle = mutated_ts.sequence_length / 2

min_dist = float('inf')
middle_variant = None
middle_pos = None

for v in mutated_ts.variants():
    dist = abs(v.site.position - genome_middle)
    if dist < min_dist:
        min_dist = dist
        middle_variant = v
        middle_pos = v.site.position

variant = middle_variant  # Use this in the rest of your code

print(f"Selected site at position: {middle_pos}")


data = []
for ind_id in sampled_inds:
    ind = mutated_ts.individual(ind_id)
    pos = ind.location[:2]
    genotype = variant.genotypes[ind.nodes[0]]
    data.append([pos[0], pos[1], genotype])

df = pd.DataFrame(data, columns=["x", "y", "genotype"])

# === Plot 1: All Genotypes Together ===
plt.figure(figsize=(8, 8))
scatter = plt.scatter(df["x"], df["y"], c=df["genotype"], cmap="tab10", s=20)
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.title("All Genotypes (Single Variant)")
plt.colorbar(scatter, label="Genotype")
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig("sampled_1000_all_genotypes.png", dpi=300)
print("Plot saved as sampled_1000_all_genotypes.png")
#plt.show()

# === Plot 2: Separate Plot For Each Genotype ===
unique_genotypes = df["genotype"].unique()

for g in unique_genotypes:
    df_g = df[df["genotype"] == g]
    plt.figure(figsize=(8, 8))
    plt.scatter(df_g["x"], df_g["y"], color="black", s=20)
    plt.xlabel("X Position")
    plt.ylabel("Y Position")
    plt.title(f"Individuals with Genotype {g}")
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    filename = f"sampled_1000_genotype_{g}.png"
    plt.savefig(filename, dpi=300)
    print(f"Plot saved as {filename}")
    #plt.show()


# # import pyslim
# # import msprime
# # import numpy as np
# # import pandas as pd
# # import tskit
# # import matplotlib.pyplot as plt

# # # === Step 1: Load .trees and Overlay Mutations ===
# # #ts = pyslim.load("output_100000.trees")




# import tskit
# import pyslim
# import msprime
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

# # === Step 1: Load .trees from SLiM ===
# ts = tskit.load("slim/out/neutral_1500.trees")

# print("Loaded tree sequence with", ts.num_individuals, "individuals")

# # === Step 2: Overlay Neutral Mutations ===
# mutation_rate = 1e-7
# #mutated_ts = msprime.mutate(ts, rate=mutation_rate, random_seed=42, keep=True)

# next_id = pyslim.next_slim_mutation_id(ts)
# mutated_ts = msprime.sim_mutations(ts,rate=mutation_rate,model=msprime.SLiMMutationModel(type=0, next_id=next_id),keep=True)


# print("Mutations overlaid with rate =", mutation_rate)

# # Optional: save mutated .trees
# mutated_ts.dump("output_100000_mutated.trees")

# # === Step 3: Sample 1000 Individuals Uniformly ===
# sample_size = 1000
# all_inds = [ind.id for ind in mutated_ts.individuals() if len(ind.nodes) > 0]
# n_sample = min(sample_size, len(all_inds))

# sampled_inds = np.random.choice(all_inds, n_sample, replace=False)

# print("Sampled", n_sample, "individuals.")

# # === Step 4: Extract Spatial Positions + Genotype ===
# variant = next(mutated_ts.variants())  # Pick first segregating site

# data = []
# for ind_id in sampled_inds:
#     ind = mutated_ts.individual(ind_id)
#     pos = ind.location[:2]  # x, y coordinates
#     genotype = variant.genotypes[ind.nodes[0]]  # Haploid genotype for color
#     data.append([pos[0], pos[1], genotype])

# df = pd.DataFrame(data, columns=["x", "y", "genotype"])

# # === Step 5: Plot ===
# plt.figure(figsize=(8, 8))
# scatter = plt.scatter(
#     df["x"], df["y"],
#     c=df["genotype"], cmap="tab10", s=20
# )

# plt.xlabel("X Position")
# plt.ylabel("Y Position")
# plt.title("Randomly Sampled 1000 Individuals (Colored by Genotype)")
# plt.colorbar(scatter, label="Genotype")
# plt.gca().set_aspect('equal')
# plt.tight_layout()

# # === Step 6: Save Plot to PNG ===
# output_plot = "sampled_1000_individuals.png"
# plt.savefig(output_plot, dpi=300)
# print("Plot saved as", output_plot)

# plt.show()
