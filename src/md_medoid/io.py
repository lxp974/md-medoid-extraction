# src/md_medoid/io.py
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, Tuple, Optional, Dict, List
import logging

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.coordinates.XTC import XTCWriter

log = logging.getLogger(__name__)


# ---------- Loading topologies and trajectories --------

# def load_universe(topology: Path | str, trajectories: Sequence[Path | str]) -> mda.Universe:
#     """
#     Create a Universe from one topology and one or more trajectories.

#     Parameters
#     ----------
#     topology : Path | str
#         Topology file (PDB/PSF/GRO/etc.)
#     trajectories : Sequence[Path | str]
#         One or more trajectory files (XTC/DCD/TRR/...). Order matters.

#     Returns
#     -------
#     MDAnalysis.Universe
#     """
#     topo = str(Path(topology))
#     trajs = [str(Path(t)) for t in trajectories]
#     if not Path(topo).exists():
#         raise FileNotFoundError(f"Topology not found: {topology}")
#     for t in trajs:
#         if not Path(t).exists():
#             raise FileNotFoundError(f"Trajectory not found: {t}")
#     log.info("Loading Universe: %s + %d traj(s)", topo, len(trajs))
#     return mda.Universe(topo, trajs)


# def load_many(pairs: Sequence[Tuple[Path | str, Path | str | Sequence[Path | str]]]) -> List[mda.Universe]:
#     """
#     Load multiple universes given (topology, trajectory_or_list) pairs.

#     Example
#     -------
#     pairs = [
#         (Path("r1.pdb"), Path("r1.xtc")),
#         (Path("r2.pdb"), [Path("r2a.xtc"), Path("r2b.xtc")]),
#     ]
#     universes = load_many(pairs)
#     """
#     universes = []
#     for top, traj in pairs:
#         trajs = traj if isinstance(traj, (list, tuple)) else [traj]
#         universes.append(load_universe(top, trajs))  # type: ignore[arg-type]
#     return universes


Pathish = str|Path
TrajLike = Pathish | Sequence[Pathish]
Pair = Tuple[Pathish, TrajLike]

def _as_paths(top: Pathish, trajs: TrajLike) -> tuple[str, list[str]]:
    topo = str(Path(top))
    traj_list = trajs if isinstance(trajs, (list, tuple)) else [trajs]
    return topo, [str(Path(t)) for t in traj_list]

def load_universes(pairs: Sequence[Pair]) -> List[mda.Universe]:
    """Load one or more Universes from (topology, trajectories) pairs."""
    universes: List[mda.Universe] = []
    for top, trajs in pairs:
        topo, ts = _as_paths(top, trajs)
        universes.append(mda.Universe(topo, ts))
    return universes


# ---------- Align trajectories to pdb reference frame  ----------

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


def align_to_reference(
    u: mda.Universe,
    ref_pdb: Path | str,
    selection: str = "protein and name CA",
    *,
    in_memory: bool = True,
    weights: str | None = None,
) -> None:
    """
    Rigidly align `u` to a reference structure using a selection (same selection on both).

    This mutates `u`'s coordinates in-place (frames are transformed). Use `in_memory=True`
    to keep results in RAM, or False to stream (slower, less memory).

    Parameters
    ----------
    u : Universe
        Mobile universe to align (trajectory loaded).
    ref_pdb : Path | str
        Reference structure (PDB) to align to.
    selection : str
        MDAnalysis selection for alignment (should exist in both).
    in_memory : bool
        If True, write transformed coords back to the trajectory in memory.
    weights : str | None
        Optional MDAnalysis weight keyword for alignment (e.g., "mass").
    """
    ref = mda.Universe(str(Path(ref_pdb)))
    _validate_selection_compatible(u, ref, selection)
    log.info("Aligning to reference %s with selection: %s", ref_pdb, selection)
    job = align.AlignTraj(u, ref, select=selection, in_memory=in_memory, weights=weights)
    job.run()


# ---------- Writing ----------

def write_aligned(
    u: mda.Universe,
    *,
    out_xtc: Path | None = None,
    out_pdb: Path | None = None,
    start: int | None = None,
    stop: int | None = None,
    stride: int = 1,
    selection: str | None = None,
    pdb_multiframe: bool = False,
) -> Dict[str, Path]:
    """
    Write aligned coordinates to disk.

    Parameters
    ----------
    u : Universe
        (Typically just aligned) Universe to write.
    out_xtc : Path | None
        Path to XTC to write (trajectory). If None, skip.
    out_pdb : Path | None
        Path to PDB to write (structure). If `pdb_multiframe=True`, writes multiple MODELs.
        If False, writes a single-frame snapshot (first selected frame).
    start, stop, stride : int | None, int | None, int
        Frame slicing, like Python's slice.
    selection : str | None
        If provided, write only this atom selection; otherwise write all atoms.
        (For portability, many workflows write *all* atoms.)
    pdb_multiframe : bool
        If True and `out_pdb` is set, write all frames to the PDB as a multi-model file.

    Returns
    -------
    dict
        Keys present for the files that were written, e.g. {"xtc": Path(...), "pdb": Path(...)}
    """
    written: Dict[str, Path] = {}
    slice_obj = slice(start, stop, stride)
    atoms = u.select_atoms(selection) if selection else u.atoms

    # Write XTC (trajectory)
    if out_xtc is not None:
        out_xtc = Path(out_xtc)
        out_xtc.parent.mkdir(parents=True, exist_ok=True)
        log.info("Writing XTC: %s", out_xtc)
        with XTCWriter(str(out_xtc), n_atoms=atoms.n_atoms) as W:
            for _ts in u.trajectory[slice_obj]:
                W.write(atoms)
        written["xtc"] = out_xtc

    # Write PDB (structure or multi-model)
    if out_pdb is not None:
        out_pdb = Path(out_pdb)
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        log.info("Writing PDB (%s): %s", "multiframe" if pdb_multiframe else "single-frame", out_pdb)
        if pdb_multiframe:
            with PDBWriter(str(out_pdb), multiframe=True) as W:
                for _ts in u.trajectory[slice_obj]:
                    W.write(atoms)
        else:
            # Write first frame of the selection/slice
            # Ensure we're at the first frame in the requested slice
            # (MDAnalysis doesn't support direct seek to slice start; iterate once)
            first = True
            for _ts in u.trajectory[slice_obj]:
                if first:
                    with PDBWriter(str(out_pdb), multiframe=False) as W:
                        W.write(atoms)
                    first = False
                    break
            if first:
                raise ValueError("No frames to write; check start/stop/stride.")
        written["pdb"] = out_pdb

    return written


# ---------- Batch helper ----------

def align_and_write_batch(
    pairs: Sequence[Tuple[Path | str, Path | str | Sequence[Path | str]]],
    ref_pdb: Path | str,
    selection: str = "protein and name CA",
    *,
    outdir: Path = Path("."),
    prefix: str = "r",
    start: int | None = None,
    stop: int | None = None,
    stride: int = 1,
    write_pdb: bool = True,
    write_xtc: bool = True,
    write_selection: str | None = None,
    pdb_multiframe: bool = False,
) -> List[Dict[str, Path]]:
    """
    For each (top, trajs) pair:
      1) load Universe
      2) align to `ref_pdb` using `selection`
      3) write aligned outputs as {prefix}{i}_aligned.(pdb|xtc)

    Returns a list (one per pair) of dicts with written paths.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written_all: List[Dict[str, Path]] = []

    for i, (top, trajs) in enumerate(pairs, start=1):
        traj_list = trajs if isinstance(trajs, (list, tuple)) else [trajs]
        u = load_universe(top, traj_list)
        align_to_reference(u, ref_pdb, selection=selection, in_memory=True)

        out_map: Dict[str, Path] = {}
        if write_xtc:
            out_map["xtc"] = outdir / f"{prefix}{i}_aligned.xtc"
        if write_pdb:
            out_map["pdb"] = outdir / f"{prefix}{i}_aligned.pdb"

        written = write_aligned(
            u,
            out_xtc=out_map.get("xtc"),
            out_pdb=out_map.get("pdb"),
            start=start,
            stop=stop,
            stride=stride,
            selection=write_selection,  # usually None → write all atoms
            pdb_multiframe=pdb_multiframe,
        )
        written_all.append(written)

    return written_all
