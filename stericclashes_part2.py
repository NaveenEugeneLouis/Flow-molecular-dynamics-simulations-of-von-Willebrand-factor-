import numpy as np

# Load clash data (either pre-summed total per frame or shape = n_frames x 5)
data = np.load("clashes_region.npy")  # or "clashes.npy"

# If the data has multiple columns (domains), sum across columns to get total clashes per frame
if data.ndim > 1 and data.shape[1] > 1:
    total_clashes_per_frame = np.sum(data, axis=1)
else:
    total_clashes_per_frame = data

# Print some basic info
print(f"Loaded clash data with shape: {data.shape}")
print("\nRaw total clashes per frame (first 10 frames):")
print(total_clashes_per_frame[:10])

# Overall statistics
print("\n--- Overall clash statistics ---")
print(f"Mean total clashes per frame = {np.mean(total_clashes_per_frame):.2f}")
print(f"Std deviation of total clashes per frame = {np.std(total_clashes_per_frame):.2f}")
print(f"Min clashes = {np.min(total_clashes_per_frame)}, Max clashes = {np.max(total_clashes_per_frame)}")
