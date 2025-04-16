import pyslim
import tskit
import numpy as np
import matplotlib.pyplot as plt

# Load tree sequence
ts = pyslim.load("output_with_mut.trees")

# For each individual, get their spatial position and genotype at a given site
positions = []
alleles = []

# Choose a focal site to color by (e.g., first segregating site)
variant = next(ts.variants())

for ind in ts.individuals():
    if len(ind.nodes) == 0:
        continue
    pos = ind.location[:2]  # x,y coordinates
    genotype = variant.genotypes[ind.nodes[0]]  # diploid: could take both nodes
    positions.append(pos)
    alleles.append(genotype)

positions = np.array(positions)

# Plotting
plt.figure(figsize=(8, 8))
scatter = plt.scatter(positions[:, 0], positions[:, 1], c=alleles, cmap="tab10", s=10)
plt.xlabel("X")
plt.ylabel("Y")
plt.colorbar(scatter, label="Allele")
plt.title("Spatial Distribution of Alleles")
plt.show()
