import pandas as pd
import numpy as np

# Load your interaction data
df = pd.read_csv("all_interactions_by_time.tsv", sep="\t")

# Define bidirectional interaction conditions
interactions = {
    "HBOND_NAIM_A1": (
        (df["type"] == "HBOND") & 
        (((df["group1"] == "NAIM") & (df["group2"] == "A1")) |
         ((df["group1"] == "A1") & (df["group2"] == "NAIM")))
    ),
    "HBOND_CAIM_A1": (
        (df["type"] == "HBOND") & 
        (((df["group1"] == "CAIM") & (df["group2"] == "A1")) |
         ((df["group1"] == "A1") & (df["group2"] == "CAIM")))
    ),
    "HBOND_NAIM_CAIM": (
        (df["type"] == "HBOND") & 
        (((df["group1"] == "NAIM") & (df["group2"] == "CAIM")) |
         ((df["group1"] == "CAIM") & (df["group2"] == "NAIM")))
    ),
    "SALTBRIDGE_NAIM_A1": (
        (df["type"] == "SALT_BRIDGE") & 
        (((df["group1"] == "NAIM") & (df["group2"] == "A1")) |
         ((df["group1"] == "A1") & (df["group2"] == "NAIM")))
    ),
    "SALTBRIDGE_CAIM_A1": (
        (df["type"] == "SALT_BRIDGE") & 
        (((df["group1"] == "CAIM") & (df["group2"] == "A1")) |
         ((df["group1"] == "A1") & (df["group2"] == "CAIM")))
    ),
    "SALTBRIDGE_NAIM_CAIM": (
        (df["type"] == "SALT_BRIDGE") & 
        (((df["group1"] == "NAIM") & (df["group2"] == "CAIM")) |
         ((df["group1"] == "CAIM") & (df["group2"] == "NAIM")))
    ),
}

# Save each interaction as .npy file
for label, condition in interactions.items():
    interaction_df = df[condition].copy()
    all_frames = range(df["frame"].min(), df["frame"].max() + 1)
    counts = interaction_df.groupby("frame").size().reindex(all_frames, fill_value=0).values
    np.save(f"{label}.npy", counts)
    print(f"✅ Saved {label}.npy with shape {counts.shape}")
