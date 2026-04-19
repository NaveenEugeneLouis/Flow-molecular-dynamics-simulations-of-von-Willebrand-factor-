# Flow-molecular-dynamics-simulations-of-von-Willebrand-factor-
Flow Molecular Dynamics (MD) Simulation Analysis Scripts

README – Molecular Dynamics Simulations of the VWF Mechanomodule Under Flow, Force, and Glycosylation States
Overview

This dataset contains all-atom molecular dynamics (MD) simulations of the von Willebrand factor (VWF) mechanomodule, designed to investigate how hydrodynamic flow, mechanical force, and glycosylation modulate its conformational dynamics and mechanostability.

Simulations were performed using GROMACS and span equilibrium, flow-driven, and force-probe regimes. The dataset enables systematic comparison between glycosylated and non-glycosylated states, with replicated simulations to ensure statistical robustness.

Simulation Details
MD engine: GROMACS v2021.2
Force field: CHARMM force field (GROMACS implementation)
Water model: TIP3P
Thermodynamic conditions (temperature and pressure): defined in the simulation protocols described in the associated publication

All simulations include periodic boundary conditions and standard molecular dynamics integration schemes.

Biological Significance

von Willebrand factor is a key mechanosensitive protein governing hemostasis. Its mechanomodule undergoes force- and flow-induced conformational transitions that regulate platelet adhesion and activation.

This dataset enables quantitative investigation of how mechanical stress and glycosylation jointly regulate structural stability and functional transitions.

Dataset Structure
1. Equilibrium MD (Free MD)

Baseline simulations used to characterize intrinsic conformational dynamics.

Simulation length: 100–500 ns
Conditions:
Glycosylated
Non-glycosylated
2. Flow MD

Simulations designed to approximate hydrodynamic flow conditions.

Parameters varied:
Initial molecular orientation
Flow velocity
Restraint configurations
Conditions:
Glycosylated
Non-glycosylated
Replicates: Each condition was simulated in triplicate
3. Steered Molecular Dynamics (SMD)

Force-probe simulations used to investigate mechanical response pathways.

External forces applied to induce conformational transitions
Conditions:
Glycosylated
Non-glycosylated
Replicates: Each condition was simulated in triplicate
File Contents
.xtc — trajectory files containing atomic coordinates over time
.pdb — reference structure of the VWF mechanomodule
topol.top — system topology defining molecular composition and parameters
Usage Example (GROMACS)
gmx trjconv -s structure.pdb -f traj.xtc -o output.pdb
Glycosylation States
With glycans: explicit glycan chains included
Without glycans: corresponding deglycosylated control systems

Notes 
No .tpr run input files are included
Parameter files (.mdp) are included and were derived from established literature protocols.
https://pubs.acs.org/doi/10.1021/acs.jpcb.4c04901
http://www.mdtutorials.com/gmx/


Suggested Applications
Conformational dynamics analysis
Glycosylation-dependent structural comparison
Flow- and force-induced transition analysis
Benchmarking of MD analysis workflows
Code Availability

Analysis scripts are available in the associated GitHub repository:
https://github.com/NaveenEugeneLouis/Flow-molecular-dynamics-simulations-of-von-Willebrand-factor-.git

Scripts include:

Intermolecular interaction analysis
Contact frequency analysis
Steric clash analysis (reference structure: PDB 1SQ0)
Citation

If you use this dataset, please cite:

Publication: https://doi.org/10.64898/2026.04.04.716521
Zenodo dataset: https://doi.org/10.5281/zenodo.19656129
Link to GitHub : https://github.com/NaveenEugeneLouis/Flow-molecular-dynamics-simulations-of-von-Willebrand-factor-.git
