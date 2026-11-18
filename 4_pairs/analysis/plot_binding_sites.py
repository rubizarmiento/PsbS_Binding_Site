import yaml
import matplotlib.pyplot as plt

# Load chain -> color mapping
with open(f"{wdir}/definitions_yaml/chain_colors.yaml", "r") as f:
    dict_chain_colors = yaml.safe_load(f)

# Chains to highlight: CP24, CP26, CP29, CP43 (both monomers)
chains = ["4", "s", "r", "c", "8", "S", "R", "C"]

fig, ax = plt.subplots(figsize=(10, 10))

# PSII background
ax.scatter(pos0[:, 0], pos0[:, 1], color="lightgrey", s=70)

# Highlighted chains
for chain in chains:
    positions_chain = dict_daft[f"chain_{chain}_seed{seeds[0]}_0001"].positions_chain
    ax.scatter(positions_chain[:, 0], positions_chain[:, 1],
               color=dict_chain_colors[chain], s=100)

# Clean layout: equal aspect, no ticks, no borders
ax.set_aspect("equal")
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

plt.savefig(f"{wdir}/4_pairs/analysis/figures/binding_sites/binding_sites.png",
            transparent=True, bbox_inches="tight", dpi=600)
plt.show()
