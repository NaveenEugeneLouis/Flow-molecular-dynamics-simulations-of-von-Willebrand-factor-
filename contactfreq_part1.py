import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import distance_array, calc_angles

# === User Input ===
topology = "TRAJECTORY.pdb"
trajectory = "TRAJECTORY.xtc"

# Load Universe
u = mda.Universe(topology, trajectory)

# Define selections (renamed)
NAIM = u.select_atoms("resid 1238-1271 and protein")
A1   = u.select_atoms("resid 1272-1458 and protein")
CAIM = u.select_atoms("resid 1459-1493 and protein")

# Donor and acceptor atom types
donor_atom_names = ["N", "NE", "ND1", "NE2", "NH1", "NH2", "OG", "OG1", "OH"]
acceptor_atom_names = ["O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"]

# === Helper: Salt bridge atoms ===
def get_salt_atoms(group):
    basic = group.select_atoms("(resname ARG and name NE NH1 NH2) or (resname LYS and name NZ)")
    acidic = group.select_atoms("(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)")
    return basic, acidic

# === Hydrogen bond detection ===
def find_hbonds(group1, group2, dist_cutoff=3.5, angle_cutoff=150):
    results = []
    for ts in u.trajectory:
        donors = group1.select_atoms("name " + " ".join(donor_atom_names))
        acceptors = group2.select_atoms("name " + " ".join(acceptor_atom_names))

        for donor in donors:
            dists = distance_array(donor.position[np.newaxis, :], acceptors.positions)[0]
            within_cutoff = np.where(dists <= dist_cutoff)[0]
            for idx in within_cutoff:
                acceptor = acceptors[idx]
                angle_rad = calc_angles(
                    donor.position[np.newaxis, :],
                    donor.position[np.newaxis, :],
                    acceptor.position[np.newaxis, :]
                )[0]
                angle_deg = np.degrees(angle_rad)
                if angle_deg >= angle_cutoff:
                    results.append([
                        ts.frame,
                        donor.resid, donor.resname, donor.name,
                        acceptor.resid, acceptor.resname, acceptor.name,
                        round(dists[idx], 2)
                    ])
    columns = ["frame", "resid1", "resname1", "atom1", "resid2", "resname2", "atom2", "distance"]
    return pd.DataFrame(results, columns=columns)

# === Salt bridge detection ===
def find_salt_bridges(group1, group2, cutoff=3.5):
    results = []
    for ts in u.trajectory:
        basic1, acidic1 = get_salt_atoms(group1)
        basic2, acidic2 = get_salt_atoms(group2)

        d1 = distance_array(basic1.positions, acidic2.positions)
        d2 = distance_array(basic2.positions, acidic1.positions)

        for i, j in zip(*np.where(d1 <= cutoff)):
            results.append([
                ts.frame,
                basic1[i].resid, basic1[i].resname, basic1[i].name,
                acidic2[j].resid, acidic2[j].resname, acidic2[j].name,
                round(d1[i, j], 2)
            ])

        for i, j in zip(*np.where(d2 <= cutoff)):
            results.append([
                ts.frame,
                basic2[i].resid, basic2[i].resname, basic2[i].name,
                acidic1[j].resid, acidic1[j].resname, acidic1[j].name,
                round(d2[i, j], 2)
            ])
    columns = ["frame", "resid1", "resname1", "atom1", "resid2", "resname2", "atom2", "distance"]
    return pd.DataFrame(results, columns=columns)

# === Run analyses ===
print("🔍 Salt bridges (NAIM ↔ A1)")
sb_NAIM_A1_df = find_salt_bridges(NAIM, A1)
sb_NAIM_A1_df.to_csv("salt_bridges_NAIM_A1.tsv", sep="\t", index=False)

print("🔍 Salt bridges (CAIM ↔ A1)")
sb_CAIM_A1_df = find_salt_bridges(CAIM, A1)
sb_CAIM_A1_df.to_csv("salt_bridges_CAIM_A1.tsv", sep="\t", index=False)

print("🔍 Salt bridges (NAIM ↔ CAIM)")
sb_NAIM_CAIM_df = find_salt_bridges(NAIM, CAIM)
sb_NAIM_CAIM_df.to_csv("salt_bridges_NAIM_CAIM.tsv", sep="\t", index=False)

print("🔍 Hydrogen bonds (NAIM ↔ A1)")
hb_NAIM_A1_df = find_hbonds(NAIM, A1)
hb_NAIM_A1_df.to_csv("hbonds_NAIM_A1.tsv", sep="\t", index=False)

print("🔍 Hydrogen bonds (CAIM ↔ A1)")
hb_CAIM_A1_df = find_hbonds(CAIM, A1)
hb_CAIM_A1_df.to_csv("hbonds_CAIM_A1.tsv", sep="\t", index=False)

print("🔍 Hydrogen bonds (NAIM ↔ CAIM)")
hb_NAIM_CAIM_df = find_hbonds(NAIM, CAIM)
hb_NAIM_CAIM_df.to_csv("hbonds_NAIM_CAIM.tsv", sep="\t", index=False)

print("✅ All interaction files saved.")
