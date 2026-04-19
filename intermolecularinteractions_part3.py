import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
from scipy.signal import savgol_filter

# === Input groups ===
input_groups = {
    "NAIM–A1": {
        "GFLOWiles": ["RUN1_SALTBRIDGE_NAIM_A1.npy", "RUN2_SALTBRIDGE_NAIM_A1.npy", "RUN3_SALTBRIDGE_NAIM_A1.npy"],
        "color": "red"
    },
    "CAIM–A1": {
        "GFLOWiles": ["RUN1_SALTBRIDGE_CAIM_A1.npy", "RUN2_SALTBRIDGE_CAIM_A1.npy", "RUN3_SALTBRIDGE_CAIM_A1.npy"],
        "color": "mediumblue"
    },
    "NAIM–CAIM": {
        "GFLOWiles": ["RUN1_SALTBRIDGE_NAIM_CAIM.npy", "RUN2_SALTBRIDGE_NAIM_CAIM.npy", "RUN3_SALTBRIDGE_NAIM_CAIM.npy"],
        "color": "orange"
    }
}

plt.figure(figsize=(12, 6))
ymax_all = 0
window_length = 15  # Must be odd
polyorder = 3

for label, info in input_groups.items():
    input_files = info["GFLOWiles"]
    color = info["color"]

    # Load and align data
    data = [np.load(f) for f in input_files]
    min_len = min(len(arr) for arr in data)
    data = [arr[:min_len] for arr in data]

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
    plt.plot(frames, smooth_mean, label=label, color=color, linewidth=2)
    plt.fill_between(frames, smooth_mean - smooth_std, smooth_mean + smooth_std,
                     color=color, alpha=0.3)

    ymax_all = max(ymax_all, np.max(smooth_mean + smooth_std))

# === Set labels and style ===
plt.xlabel("GFLOWrame")
plt.ylabel("Number of Hydrogen Bonds")
plt.title("Hydrogen Bonds per Frame (Smoothed with Savitzky–Golay Filter)")
plt.legend()

# === Remove grid lines ===
plt.grid(False)

# === Force y-axis ticks to be whole numbers ===
ymax_int = math.ceil(ymax_all)
plt.ylim(0, ymax_int)
plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(1))  # Every 1 unit

plt.tight_layout()
plt.savefig("SALTBRIDGE_curly_smoothed_plot.png", dpi=300)
plt.show()
