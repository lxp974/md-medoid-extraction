#!/usr/bin/env python3

"""
Test script to wire together the modules and to check if each
module produces the desired outputs

"""

# VERY minimal, single pair only (no validation)
import argparse
from pathlib import Path
import MDAnalysis as mda
from md_medoid.io import align_to_reference, write_aligned

p = argparse.ArgumentParser()
p.add_argument("-t", "--topology", required=True)
p.add_argument("-x", "--traj", required=True)
p.add_argument("--ref", required=True)
p.add_argument("--sel-align", default="protein and name CA")
p.add_argument("-o", "--outdir", default="aligned")
args = p.parse_args()

u = mda.Universe(args.topology, [args.traj])
align_to_reference(u, args.ref, selection=args.sel_align)
outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
write_aligned(u, out_xtc=outdir / "aligned.xtc", out_pdb=outdir / "aligned.pdb")
print("done")
