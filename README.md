# md-medoid-extraction
**This repo is a work in progress - actively converting a research notebook into a Python package.**

Tools to cluster MD simulation trajectories and extract representative medoid frames for virtual screening.
- Takes MD simulations with n reps and aligns to a reference structure.
- Calculates RMSD of binding pocket residues
- PCA on binding pocket selected atoms
- Clusters using HDBSCAN in PCA space
- Extracts representative medoid frames from clusters 

