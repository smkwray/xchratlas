from __future__ import annotations

from dataclasses import dataclass

from .constants import PAR1_END, PAR1_START, PAR2_END, PAR2_START


@dataclass(frozen=True)
class XRegion:
    name: str
    start: int
    end: int


PAR1 = XRegion("PAR1", PAR1_START, PAR1_END)
NONPAR = XRegion("nonPAR", PAR1_END + 1, PAR2_START - 1)
PAR2 = XRegion("PAR2", PAR2_START, PAR2_END)


def classify_x_region(pos: int) -> str:
    if PAR1.start <= pos <= PAR1.end:
        return PAR1.name
    if PAR2.start <= pos <= PAR2.end:
        return PAR2.name
    if NONPAR.start <= pos <= NONPAR.end:
        return NONPAR.name
    return "outside_chrX_modeled_range"


def classify_interval_region(start: int, end: int) -> str:
    if end < start:
        start, end = end, start
    regions = {classify_x_region(start), classify_x_region(end)}
    if len(regions) == 1:
        return next(iter(regions))
    return "boundary_crossing"


def within_same_x_region(pos_a: int, pos_b: int) -> bool:
    return classify_x_region(pos_a) == classify_x_region(pos_b)
