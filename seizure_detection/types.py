from dataclasses import dataclass
from typing import Optional


@dataclass
class SeizureEvent:
    start_sec: float
    end_sec: float
    duration_sec: float
    gamma_peak: float
    pis_found: bool
    pis_start_sec: Optional[float]
    pis_end_sec: Optional[float]


@dataclass
class PeakMaskSeg:
    start_idx: int
    end_idx: int
    start_sec: float
    end_sec: float
    duration_sec: float
    peak_abs: float
