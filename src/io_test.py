# src/md_medoid/io.py
from __future__ import annotations
from pathlib import Path
from typing import Sequence, Iterable, Tuple, Union, List
import logging
import MDAnalysis as mda

log = logging.getLogger(__name__)

Pathish = Union[str, Path]
TrajLike = Union[Pathish, Sequence[Pathish]]
Pair = Tuple[Pathish, TrajLike]


def _as_str_paths(top: Pathish, trajs: TrajLike) -> tuple[str, list[str]]:
    topo = str(Path(top))
    traj_list = trajs if isinstance(trajs, (list, tuple)) else [trajs]
    trajs_str = [str(Path(t)) for t in traj_list]
    if not Path(topo).exists():
        raise FileNotFoundError(f"Topology not found: {top}")
    for t in trajs_str:
        if not Path(t).exists():
            raise FileNotFoundError(f"Trajectory not found: {t}")
    return topo, trajs_str


def load_pair(topology: Pathish, trajectories: TrajLike) -> mda.Universe:
    """Load a single Universe from one topology and one-or-more trajectory files."""
    topo, trajs = _as_str_paths(topology, trajectories)
    log.info("Loading Universe: %s + %d traj(s)", topo, len(trajs))
    return mda.Universe(topo, trajs)


def load_pairs(pairs: Sequence[Pair]) -> List[mda.Universe]:
    """
    Load multiple Universes from (topology, trajectory/trajectories) pairs.

    Examples
    --------
    # 1 xtc per run
    ulist = load_pairs([
        ("Run1/md_prot.gro", "Run1/md_prot_cent_ts100.xtc"),
        ("Run2/md_prot.gro", "Run2/md_prot_cent_ts100.xtc"),
    ])

    # multiple xtcs per run (concatenated)
    ulist = load_pairs([
        ("Run1/md_prot.gro", ["Run1/part1.xtc", "Run1/part2.xtc"]),
        ("Run2/md_prot.gro", ["Run2/part1.xtc", "Run2/part2.xtc"]),
    ])
    """
    universes: List[mda.Universe] = []
    for top, trajs in pairs:
        universes.append(load_pair(top, trajs))
    return universes
