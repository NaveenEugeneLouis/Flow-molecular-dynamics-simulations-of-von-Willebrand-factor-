import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === USER CONFIGURATION ===
input_files = ["RUN1.tsv", "RUN2.tsv", "RUN3.tsv"]
output_image = "contact_heatmap.png"
cutoff_frequency = 0.10

# === Charged Residue Color Map ===
charge_colors = {
    "ASP": "red", "GLU": "red",
    "ARG": "blue", "LYS": "blue", "HIP": "blue"
}

# === Load and combine data ===
dfs = []
total_frames = 0
for file in input_files:
    df = pd.read_csv(file, sep="\t", skiprows=1, header=None, 
                     names=["frame", "resid1", "resname1", "atom1", "resid2", "resname2", "atom2", "distance"],
                     dtype={"frame": int, "resid1": int, "resname1": str, "atom1": str,
                            "resid2": int, "resname2": str, "atom2": str, "distance": float})
    max_frame = df["frame"].max()
    df["frame"] += total_frames  # Re-index frames to make them continuous across files
    total_frames += max_frame + 1
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

# === Processing ===
df["label1"] = df["resname1"] + df["resid1"].astype(str)  # B1
df["label2"] = df["resname2"] + df["resid2"].astype(str)  # A2

# Compute frequency of contacts
contact_counts = df.groupby(["label1", "label2"]).size().reset_index(name="count")
contact_counts["frequency"] = contact_counts["count"] / total_frames
filtered = contact_counts[contact_counts["frequency"] >= cutoff_frequency]

# Build frequency matrix
x_labels = sorted(filtered["label2"].unique(), key=lambda x: int(''.join(filter(str.isdigit, x))))
y_labels = sorted(filtered["label1"].unique(), key=lambda x: int(''.join(filter(str.isdigit, x))))
matrix = pd.DataFrame(0.0, index=y_labels, columns=x_labels)

for _, row in filtered.iterrows():
    matrix.loc[row["label1"], row["label2"]] = row["frequency"]

# === Plotting ===
plt.figure(figsize=(12, 10))
ax = sns.heatmap(
    matrix,
    cmap="hot",
    linewidths=0.5,
    linecolor="gray",
    cbar_kws={'label': 'Contact Frequency'},
    vmin=cutoff_frequency,
    vmax=1.0,
    square=True
)

# Set axis labels
ax.set_xticks(np.arange(len(x_labels)) + 0.5)
ax.set_yticks(np.arange(len(y_labels)) + 0.5)
ax.set_xticklabels(x_labels, rotation=90, fontname="Arial", fontsize=12)
ax.set_yticklabels(y_labels, rotation=0, fontname="Arial", fontsize=12)

# Color tick labels by charge
def get_label_color(label):
    for res in charge_colors:
        if label.startswith(res):
            return charge_colors[res]
    return "black"

for tick in ax.get_xticklabels():
    tick.set_color(get_label_color(tick.get_text()))
for tick in ax.get_yticklabels():
    tick.set_color(get_label_color(tick.get_text()))

# Customize colorbar
cbar = ax.collections[0].colorbar
cbar.set_label('Contact Frequency', fontname='Arial', fontsize=12)
cbar.ax.tick_params(labelsize=12)
for label in cbar.ax.get_yticklabels():
    label.set_fontname('Arial')

# Optional: Glycan marker
glycan_start = next((i for i, x in enumerate(x_labels) if "GalNAc" in x), None)
if glycan_start is not None:
    ax.axvline(x=glycan_start, color='black', linestyle='--', linewidth=1)

plt.title("Figure S10: Contact Frequency Heatmap (B1 ↔ A2)", fontsize=12, fontname="Arial")
plt.tight_layout()
plt.savefig(output_image, dpi=300)
plt.show()
