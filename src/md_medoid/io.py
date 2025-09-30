# src/md_medoid/io.py
from __future__ import annotations
from pathlib import Path
import logging
from typing import Iterable, Sequence, Tuple, Optional, Dict, List

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.coordinates.XTC import XTCWriter



log = logging.getLogger(__name__)

# ---------- Set paths --------
Pathish = str|Path
TrajLike = Pathish | Sequence[Pathish]
Pair = Tuple[Pathish, TrajLike]

def _as_paths(top: Pathish, trajs: TrajLike) -> tuple[str, list[str]]:
    topo = str(Path(top))
    traj_list = trajs if isinstance(trajs, (list, tuple)) else [trajs]
    return topo, [str(Path(t)) for t in traj_list]


# ---------- Loading topologies and trajectories --------
def load_universes(pairs: Sequence[Pair]) -> List[mda.Universe]:
    """Load one or more Universes from (topology, trajectories) pairs."""
    universes: List[mda.Universe] = []
    for top, trajs in pairs:
        topo, ts = _as_paths(top, trajs)
        universes.append(mda.Universe(topo, ts))
    return universes

# ---------- Align domains --------
# Check universes are populated
def _validate_selection_compatible(u: mda.Universe, ref: mda.Universe, selection: str) -> None:
    """
    Ensure the selection matches a non-zero, equal number of atoms in both mobile and reference.
    """
    mob = u.select_atoms(selection)
    rsel = ref.select_atoms(selection)
    if mob.n_atoms == 0:
        raise ValueError(f"Selection matched 0 atoms in mobile universe: {selection!r}")
    if rsel.n_atoms == 0:
        raise ValueError(f"Selection matched 0 atoms in reference: {selection!r}")
    if mob.n_atoms != rsel.n_atoms:
        raise ValueError(
            f"Selection size mismatch: mobile={mob.n_atoms}, reference={rsel.n_atoms} for {selection!r}"
        )
    
def align_to_ref(
        u: mda.Universe,
        ref_pdb: Path | str,
        selection : str | None = None,
        in_memory: bool =True,
    ):
    """
    Align universe 'u' to reference structure using selection set of atoms.
    
    Parameters
    -----
    u : Universe
        Mobile universe to be aligned.
    ref_pdb : Path | str
        Reference structure to align to (PDB)
    selection: str
        Selection of atoms for alignment.
    in_memory : bool
        If True, write transformed coordinates back to the trajectory in memory.

    """

    # Load reference 
    ref = mda.Universe(str(Path(ref_pdb)))
    _validate_selection_compatible(u, ref, selection)
    log.info("Aligning trajectory to reference %s with selection %s", ref_pdb, selection)
    align = align.AlignTraj(u, ref, select=selection, in_memory=in_memory)
    align.run()

    # ---------- Write aligned trajectories --------

    def write_aligned(
            u: mda.Universe,
            out_xtc: Path | None = None,
            out_pdb: Path | None = None,
            start: int | None = None,
            stop: int | None = None,
            stride: int = 1,
            selection: str | None = None,
            pdb_multiframe: bool = False,
        ) -> Dict[str, Path]:
        """
        Write the aligned coordinates/trajectories to disk
        
        Parameters
        ------
        u : Universe
            Mobile universe to be aligned.   
        out_xtc: 
            Path to xtc output to write 
        out_pdb: 
            Path to pdb output to write
        start: int
            Frame to start writing from
        stop: int
            Frame to stop writing
        stride: int 
            Write ever 'stride' frames
        selection: str
            Selection of atoms to write
        pdb_multiframe: bool
            If True, and 'out_pdb' is set, write all frames to the PDB uas a multi-model file.
        """
        # Dictionary of written files
        written: Dict[str, Path] ={}


        # Write XTC
        if out_xtc is not None:
            out_xtc=Path(out_xtc)
            out_xtc.parent.mkdir(parents=True, exist_ok=True)
            log.info("Writing xtc file: %s", out_xtc)
            with XTCWriter(str(out_xtc), n_atoms=u.atoms.n_atoms) as w:
                for ts in u.trajectory[start:stop:stride]:
                    if selection is not None:
                        w.write(u.select_atoms(selection))
            written['xtc'] = out_xtc

        # Write PDB
        if out_pdb is not None:
            out_pdb = Path(out_pdb)
            out_pdb.parent.mkdir(parents=True, exist_ok=True)
            log.info("Writting pdb file: %s", out_pdb)
            u.atoms.write(out_pdb)

        else: 
            # Write first frame
            with PDBWriter(str(out_pdb), multiframe=pdb_multiframe) as w:
                for ts in u.trajectory[start:stop:stride]:
                    if selection is not None:
                        w.write(u.select_atoms(selection))
                    else:
                        w.write(u.atoms)
            written['pdb'] = out_pdb