import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import savgol_filter

# ==========================
# Input files per condition
# ==========================
conditions = {
    "DD3NFP": {
        "files": [
            "DD3NFPclashes1.npy",
            "DD3NFPclashes2.npy",
            "DD3NFPclashes3.npy",
        ],
        "color": "red",
    },
    "DD3OG": {
        "files": [
            "DD3OGclashes1.npy",
            "DD3OGclashes2.npy",
            "DD3OGclashes3.npy",
        ],
        "color": "blue",
    },
    "DD3": {
        "files": [
            "DD3clashes1.npy",
            "DD3clashes2.npy",
            "DD3clashes3.npy",
        ],
        "color": "green",
    },
    "ProVWF": {
        "files": [
            "ProVWFclashes1.npy",
            "ProVWFclashes2.npy",
            "ProVWFclashes3.npy",
        ],
        "color": "purple",
    },
}

# ==========================
# Savitzky–Golay parameters
# ==========================
window_length = 15  # must be odd
polyorder = 3

# ==========================
# Load ALL data first (for global normalization)
# ==========================
all_data = []

for cond in conditions.values():
    for f in cond["files"]:
        arr = np.load(f)
        if arr.ndim > 1 and arr.shape[1] > 1:
            arr = np.sum(arr, axis=1)
        all_data.append(arr)

# Align lengths globally
min_len = min(len(arr) for arr in all_data)
all_data = [arr[:min_len] for arr in all_data]

# Global min/max for 0–100% normalization
global_min = min(arr.min() for arr in all_data)
global_max = max(arr.max() for arr in all_data)

if global_max == global_min:
    raise ValueError("All clash values identical — cannot normalize.")

# ==========================
# Plot
# ==========================
plt.figure(figsize=(12, 6))
frames = np.arange(min_len)

for label, cond in conditions.items():
    reps = []

    for f in cond["files"]:
        arr = np.load(f)
        if arr.ndim > 1 and arr.shape[1] > 1:
            arr = np.sum(arr, axis=1)
        arr = arr[:min_len]

        # Normalize to 0–100%
        norm = (arr - global_min) / (global_max - global_min) * 100.0
        reps.append(norm)

    reps = np.vstack(reps)
    mean_vals = np.mean(reps, axis=0)
    std_vals = np.std(reps, axis=0)

    # Smooth
    if len(mean_vals) >= window_length:
        mean_vals = savgol_filter(mean_vals, window_length, polyorder)
        std_vals = savgol_filter(std_vals, window_length, polyorder)

    print(f"{label}: Mean = {np.mean(mean_vals):.2f} %, Std = {np.mean(std_vals):.2f} %")

    plt.plot(frames, mean_vals, label=label, color=cond["color"], linewidth=2)
    plt.fill_between(
        frames,
        mean_vals - std_vals,
        mean_vals + std_vals,
        color=cond["color"],
        alpha=0.25,
    )

# ==========================
# Formatting
# ==========================
plt.xlabel("Frame")
plt.ylabel("Steric clashes (%)")
plt.title("Normalized steric clashes vs frame (0–100%)")
plt.legend()
plt.grid(False)
plt.ylim(0, 100)
plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig("Steric_clashes_normalized_smoothed.png", dpi=300)
plt.show()
