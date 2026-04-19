import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import distance_array, calc_angles

# === User Input ===
topology = "Trajectory.pdb"
trajectory = "Trajectory.xtc"

# Load Universe
u = mda.Universe(topology, trajectory)

# Define selections
NAIM = u.select_atoms("resid 1238-1271 and protein")
A1 = u.select_atoms("resid 1272-1458 and protein")
CAIM = u.select_atoms("resid 1459-1493 and protein")

donor_atom_names = ["N", "NE", "ND1", "NE2", "NH1", "NH2", "OG", "OG1", "OH"]
acceptor_atom_names = ["O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"]

interaction_pairs = [
    ("NAIM", "A1", NAIM, A1),
    ("NAIM", "CAIM", NAIM, CAIM),
    ("CAIM", "A1", CAIM, A1),
]

# === Helper: Salt bridge atoms ===
def get_salt_atoms(group):
    basic = group.select_atoms("(resname ARG and name NE NH1 NH2) or (resname LYS and name NZ)")
    acidic = group.select_atoms("(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)")
    return basic, acidic

# === Store all results ===
all_results = []

# === Frame loop ===
for ts in u.trajectory:
    frame = ts.frame

    for label1, label2, group1, group2 in interaction_pairs:
        # Hydrogen bonds
        donors = group1.select_atoms("name " + " ".join(donor_atom_names))
        acceptors = group2.select_atoms("name " + " ".join(acceptor_atom_names))

        if len(donors) and len(acceptors):
            dists = distance_array(donors.positions, acceptors.positions)
            for i, j in zip(*np.where(dists <= 3.5)):
                angle_rad = calc_angles(
                    donors[i].position[np.newaxis, :],
                    donors[i].position[np.newaxis, :],
                    acceptors[j].position[np.newaxis, :]
                )[0]
                if np.degrees(angle_rad) >= 150:
                    all_results.append([
                        frame, "HBOND", label1, label2,
                        donors[i].resid, donors[i].resname, donors[i].name,
                        acceptors[j].resid, acceptors[j].resname, acceptors[j].name,
                        round(dists[i, j], 2)
                    ])

        # Salt bridges
        basic1, acidic1 = get_salt_atoms(group1)
        basic2, acidic2 = get_salt_atoms(group2)

        d1 = distance_array(basic1.positions, acidic2.positions)
        d2 = distance_array(basic2.positions, acidic1.positions)

        for i, j in zip(*np.where(d1 <= 3.5)):
            all_results.append([
                frame, "SALT_BRIDGE", label1, label2,
                basic1[i].resid, basic1[i].resname, basic1[i].name,
                acidic2[j].resid, acidic2[j].resname, acidic2[j].name,
                round(d1[i, j], 2)
            ])

        for i, j in zip(*np.where(d2 <= 3.5)):
            all_results.append([
                frame, "SALT_BRIDGE", label2, label1,
                basic2[i].resid, basic2[i].resname, basic2[i].name,
                acidic1[j].resid, acidic1[j].resname, acidic1[j].name,
                round(d2[i, j], 2)
            ])

# === Save to DataFrame and TSV ===
columns = [
    "frame", "type", "group1", "group2",
    "resid1", "resname1", "atom1",
    "resid2", "resname2", "atom2",
    "distance"
]

df = pd.DataFrame(all_results, columns=columns)
df.to_csv("all_interactions_by_time.tsv", sep="\t", index=False)

print("✅ All interactions written to all_interactions_by_time.tsv")
