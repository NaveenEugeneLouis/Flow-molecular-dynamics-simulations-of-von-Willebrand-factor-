import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import savgol_filter

# === Input files and plot color ===
input_files = ["Run1clashes.npy", "Run2clashes.npy", "Run3clashes.npy"]  # Replace with your files
color = "red"
label = "Total clashes"

# Savitzky–Golay smoothing parameters
window_length = 15  # Must be odd
polyorder = 3

# === Load and process data ===
data = []
for f in input_files:
    arr = np.load(f)  # shape = (n_frames, 5) or (n_frames,)
    if arr.ndim > 1 and arr.shape[1] > 1:
        total = np.sum(arr, axis=1)  # sum across all domains
    else:
        total = arr
    data.append(total)

# Align lengths in case of different trajectory lengths
min_len = min(len(arr) for arr in data)
data = [arr[:min_len] for arr in data]

# Stack replicates and compute mean & std
arr = np.vstack(data)
mean_vals = np.mean(arr, axis=0)
std_vals = np.std(arr, axis=0)

# Smooth the data
if len(mean_vals) >= window_length:
    smooth_mean = savgol_filter(mean_vals, window_length, polyorder)
    smooth_std = savgol_filter(std_vals, window_length, polyorder)
else:
    smooth_mean = mean_vals
    smooth_std = std_vals

print(f"{label}: Mean = {np.mean(mean_vals):.2f}, Std = {np.mean(std_vals):.2f}")

frames = np.arange(len(smooth_mean))

# Dynamic Y-axis limits
y_min = np.min(smooth_mean - smooth_std)
y_max = np.max(smooth_mean + smooth_std)

# === Plot ===
plt.figure(figsize=(12, 6))
plt.plot(frames, smooth_mean, label=label, color=color, linewidth=2)
plt.fill_between(frames, smooth_mean - smooth_std, smooth_mean + smooth_std,
                 color=color, alpha=0.3)

plt.xlabel("Frame")
plt.ylabel("Number of Clashes (Total)")
plt.title(f"Total Clashes per Frame\nSmoothed with Savitzky–Golay Filter")
plt.legend()
plt.grid(False)
plt.ylim(y_min - 10, y_max + 10)
plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
out_filename = f"Clashes_smoothed_total.png"
plt.savefig(out_filename, dpi=300)
plt.show()
