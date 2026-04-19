import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np

# === Load dynamic trajectory (A1 domain in SMD, with glycans) ===
u = mda.Universe("Trajectory.pdb", "Trajectory.xtc")

# === Load static structure (1SQ0 with GPIbα) ===
gpib_static = mda.Universe("1sq0.pdb")

# === Alignment atoms: residues 1290–1430 (heavy atoms only) ===
mov_sel_static = gpib_static.select_atoms("resid 1290-1430 and not name H*").sort()
ref_sel = u.select_atoms("resid 1290-1430 and not name H*").sort()

if mov_sel_static.n_atoms != ref_sel.n_atoms:
    raise ValueError(
        f"Mismatch in atom count for alignment: "
        f"static={mov_sel_static.n_atoms}, traj={ref_sel.n_atoms}"
    )

# === Static GPIbα residues for clash detection (1–265) ===
gpib_all_static = gpib_static.select_atoms("resid 1-265 and not name H*")

# === Glycan resnames ===
glycan_resnames = ["AGAL", "BGLC", "AFUC", "ANE5", "BGAL", "BMAN", "AMAN"]
glycan_str = " or ".join([f"resname {g}" for g in glycan_resnames])
heavy = "not name H*"

# === Single dynamic region selection
#     (resid 764–1271 + nearby glycans AND resid 1493–1873 + nearby glycans)
Region_sel = u.select_atoms(
    f"((resid 764-1271 or resid 1493-1873) "
    f"or ({glycan_str} and around 5 (resid 764-1271 or resid 1493-1873))) "
    f"and {heavy}"
)

# === Initialize output ===
all_clashes = []

# === Per-frame loop ===
for i, ts in enumerate(u.trajectory[0:3001]):
    try:
        # Alignment positions
        mov_pos = mov_sel_static.positions.copy()
        ref_pos = ref_sel.positions.copy()
        gpib_pos = gpib_all_static.positions.copy()

        # Compute rotation + translation
        R, rmsd = rotation_matrix(mov_pos, ref_pos)
        mov_com = mov_pos.mean(axis=0)
        ref_com = ref_pos.mean(axis=0)

        aligned_gpib_pos = np.dot((gpib_pos - mov_com), R.T) + ref_com

        # Compute clashes
        clash_count = np.sum(
            distance_array(Region_sel.positions, aligned_gpib_pos) < 3.0
        )
        all_clashes.append(clash_count)

        if i % 100 == 0:
            print(f"Processed frame {ts.frame}, clashes = {clash_count}")

    except Exception as e:
        print(f"⚠️ Skipping frame {ts.frame} due to error: {e}")
        continue

# === Save results ===
all_clashes = np.array(all_clashes)
np.save("clashes_region.npy", all_clashes)
print(
    f"✅ Clash analysis complete. "
    f"Saved clashes_region.npy with shape {all_clashes.shape}"
)
