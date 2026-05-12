from typing import List, Optional, Tuple

import numpy as np
from scipy import signal

from .types import PeakMaskSeg, SeizureEvent


def _find_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Return inclusive start/end index runs where mask is True."""
    mask = np.asarray(mask, dtype=bool)

    if mask.size == 0 or not np.any(mask):
        return []

    edges = np.diff(mask.astype(np.int8))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0]

    if mask[0]:
        starts = np.r_[0, starts]
    if mask[-1]:
        ends = np.r_[ends, len(mask) - 1]

    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def mask_super_peaks_adaptive(
    x: np.ndarray,
    fs: float,
    peak_thr: float = 0.0015,
    recover_thr: Optional[float] = None,
    stable_sec: float = 0.2,
    max_expand_sec: float = 2.0,
    merge_gap_sec: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, List[PeakMaskSeg]]:
    """
    Mask very large EEG peaks by expanding around threshold-crossing regions.

    Returns
    -------
    x_out:
        Signal with masked regions replaced by NaN.
    mask:
        Boolean mask. True means the sample was masked.
    segs:
        Metadata for each masked segment.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)

    if n == 0:
        return x.copy(), np.zeros(0, dtype=bool), []

    if recover_thr is None:
        recover_thr = peak_thr / 3.0

    stable_n = max(1, int(round(stable_sec * fs)))
    max_expand_n = max(1, int(round(max_expand_sec * fs)))
    merge_gap_n = max(0, int(round(merge_gap_sec * fs)))

    peak_mask = np.abs(x) > peak_thr
    runs = _find_runs(peak_mask)

    if not runs:
        return x.copy(), np.zeros(n, dtype=bool), []

    segs_idx: List[Tuple[int, int, float]] = []

    for s, e in runs:
        peak_abs = float(np.max(np.abs(x[s : e + 1])))

        left = s
        steps = 0
        while left > 0 and steps < max_expand_n:
            l2 = max(0, left - stable_n)
            if np.all(np.abs(x[l2:left]) < recover_thr):
                break
            left -= 1
            steps += 1

        right = e
        steps = 0
        while right < n - 1 and steps < max_expand_n:
            r2 = min(n, right + stable_n + 1)
            if np.all(np.abs(x[right + 1 : r2]) < recover_thr):
                break
            right += 1
            steps += 1

        segs_idx.append((left, right, peak_abs))

    segs_idx.sort(key=lambda z: z[0])

    merged: List[Tuple[int, int, float]] = []
    cur_s, cur_e, cur_pk = segs_idx[0]

    for s, e, pk in segs_idx[1:]:
        if s <= cur_e + merge_gap_n:
            cur_e = max(cur_e, e)
            cur_pk = max(cur_pk, pk)
        else:
            merged.append((cur_s, cur_e, cur_pk))
            cur_s, cur_e, cur_pk = s, e, pk

    merged.append((cur_s, cur_e, cur_pk))

    mask = np.zeros(n, dtype=bool)
    segs: List[PeakMaskSeg] = []

    for s, e, pk in merged:
        s = max(0, int(s))
        e = min(n - 1, int(e))

        if e < s:
            continue

        mask[s : e + 1] = True
        segs.append(
            PeakMaskSeg(
                start_idx=s,
                end_idx=e,
                start_sec=s / fs,
                end_sec=e / fs,
                duration_sec=(e - s + 1) / fs,
                peak_abs=pk,
            )
        )

    x_out = x.copy()
    x_out[mask] = np.nan

    return x_out, mask, segs


def gamma_power_series_stft(
    x: np.ndarray,
    fs: float,
    gamma_band: Tuple[float, float] = (20.0, 50.0),
    win_sec: float = 1.0,
    step_sec: float = 1.0,
    detrend: str = "constant",
    nan_max_frac: float = 0.2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute gamma-band power over time using sliding FFT windows.

    NaN handling:
    - If NaN fraction in one window > nan_max_frac, gamma power is NaN.
    - Otherwise, NaNs are linearly interpolated inside that window.
    """
    x = np.asarray(x, dtype=float)

    nperseg = int(round(win_sec * fs))
    hop = int(round(step_sec * fs))
    n = len(x)

    if nperseg <= 0 or hop <= 0:
        raise ValueError("win_sec and step_sec must produce positive sample counts.")

    if n < nperseg:
        return np.array([], dtype=float), np.array([], dtype=float)

    freqs = np.fft.rfftfreq(nperseg, d=1.0 / fs)
    fmask = (freqs >= gamma_band[0]) & (freqs <= gamma_band[1])

    if not np.any(fmask):
        raise ValueError("gamma_band is outside available frequency bins.")

    win = signal.windows.hann(nperseg, sym=False)

    t_list: List[float] = []
    gp_list: List[float] = []

    for start in range(0, n - nperseg + 1, hop):
        seg = x[start : start + nperseg]
        nan_mask = ~np.isfinite(seg)
        nan_frac = float(np.mean(nan_mask))

        center = (start + nperseg / 2) / fs
        t_list.append(center)

        if nan_frac > nan_max_frac:
            gp_list.append(np.nan)
            continue

        if nan_frac > 0:
            idx = np.arange(nperseg)
            good = ~nan_mask

            if np.sum(good) < 2:
                gp_list.append(np.nan)
                continue

            seg = seg.copy()
            seg[nan_mask] = np.interp(idx[nan_mask], idx[good], seg[good])

        if detrend is not None:
            seg = signal.detrend(seg, type=detrend)

        segw = seg * win
        fft_vals = np.fft.rfft(segw)
        power = np.abs(fft_vals) ** 2

        gp_list.append(float(np.sum(power[fmask])))

    return np.asarray(t_list), np.asarray(gp_list)


def find_excursions_on_series(
    series: np.ndarray,
    thr: float,
    min_len: int,
) -> List[Tuple[int, int]]:
    """
    Find contiguous segments where series > threshold.

    NaNs are treated as segment breakers.
    """
    series = np.asarray(series, dtype=float)
    above = (series > thr) & np.isfinite(series)

    if not np.any(above):
        return []

    runs = _find_runs(above)
    return [(s, e) for s, e in runs if (e - s + 1) >= min_len]


def pis_test_in_post_window(
    gamma_series: np.ndarray,
    anchor_idx: int,
    search_len: int = 90,
    min_pis_len: int = 6,
    ratio: float = 100.0,
) -> Tuple[bool, Optional[int], Optional[int]]:
    """
    Test whether a post-ictal suppression-like low-gamma period appears
    after a gamma excursion.
    """
    g = np.asarray(gamma_series, dtype=float)
    n = len(g)

    i0 = int(anchor_idx)
    i1 = min(n, i0 + int(search_len))

    if i1 - i0 < min_pis_len:
        return False, None, None

    seg = g[i0:i1].copy()
    finite = np.isfinite(seg)

    if np.sum(finite) < max(5, min_pis_len):
        return False, None, None

    vals = np.maximum(seg[finite], 1e-12)
    vals_sorted = np.sort(vals)

    k = min(3, len(vals_sorted))
    mean_min = float(np.mean(vals_sorted[:k]))
    mean_max = float(np.mean(vals_sorted[-k:]))

    if mean_max / mean_min < ratio:
        return False, None, None

    low_thr = mean_max / ratio
    low = (seg <= low_thr) & np.isfinite(seg)

    if not np.any(low):
        return False, None, None

    for s, e in _find_runs(low):
        if e - s + 1 >= min_pis_len:
            return True, i0 + int(s), i0 + int(e)

    return False, None, None


def detect_seizures_gamma_pis(
    x_full: np.ndarray,
    fs: float,
    x_baseline: np.ndarray,
    gamma_band: Tuple[float, float] = (20.0, 50.0),
    win_sec: float = 1.0,
    step_sec: float = 1.0,
    thr_sd: float = 3.0,
    min_excursion_sec: float = 5.0,
    pis_search_sec: float = 90.0,
    min_pis_sec: float = 6.0,
    pis_ratio: float = 100.0,
    enable_peak_mask: bool = True,
    peak_thr: float = 0.0015,
    recover_thr: Optional[float] = None,
    stable_sec: float = 0.2,
    max_expand_sec: float = 2.0,
    merge_gap_sec: float = 0.05,
    nan_max_frac: float = 0.2,
) -> List[SeizureEvent]:
    """
    Main seizure detector.

    Pipeline:
    1. Optionally mask super-peaks.
    2. Compute gamma power series.
    3. Use baseline gamma mean + SD threshold.
    4. Find sustained gamma excursions.
    5. Keep excursions with post-ictal suppression-like gamma drop.
    """
    x_full = np.asarray(x_full, dtype=float)
    x_baseline = np.asarray(x_baseline, dtype=float)

    if enable_peak_mask:
        x_full, _, _ = mask_super_peaks_adaptive(
            x_full,
            fs,
            peak_thr=peak_thr,
            recover_thr=recover_thr,
            stable_sec=stable_sec,
            max_expand_sec=max_expand_sec,
            merge_gap_sec=merge_gap_sec,
        )

        x_baseline, _, _ = mask_super_peaks_adaptive(
            x_baseline,
            fs,
            peak_thr=peak_thr,
            recover_thr=recover_thr,
            stable_sec=stable_sec,
            max_expand_sec=max_expand_sec,
            merge_gap_sec=merge_gap_sec,
        )

    t_full, g_full = gamma_power_series_stft(
        x_full,
        fs,
        gamma_band=gamma_band,
        win_sec=win_sec,
        step_sec=step_sec,
        nan_max_frac=nan_max_frac,
    )

    _, g_base = gamma_power_series_stft(
        x_baseline,
        fs,
        gamma_band=gamma_band,
        win_sec=win_sec,
        step_sec=step_sec,
        nan_max_frac=nan_max_frac,
    )

    g_base_finite = g_base[np.isfinite(g_base)]

    if g_base_finite.size < 5 or g_full.size == 0:
        return []

    mu = float(np.mean(g_base_finite))
    sd = float(np.std(g_base_finite, ddof=1))

    threshold = mu + thr_sd * sd

    series_hz = 1.0 / step_sec
    min_excursion_len = int(np.ceil(min_excursion_sec * series_hz))
    pis_search_len = int(np.ceil(pis_search_sec * series_hz))
    min_pis_len = int(np.ceil(min_pis_sec * series_hz))

    candidate_segments = find_excursions_on_series(
        g_full,
        thr=threshold,
        min_len=min_excursion_len,
    )

    events: List[SeizureEvent] = []

    for s_idx, e_idx in candidate_segments:
        pis_found, pis_s, pis_e = pis_test_in_post_window(
            g_full,
            anchor_idx=s_idx,
            search_len=pis_search_len,
            min_pis_len=min_pis_len,
            ratio=pis_ratio,
        )

        if not pis_found:
            continue

        start_sec = float(t_full[s_idx])
        end_sec = float(t_full[e_idx])
        duration_sec = max(0.0, end_sec - start_sec + step_sec)

        segment_gamma = g_full[s_idx : e_idx + 1]
        segment_gamma = segment_gamma[np.isfinite(segment_gamma)]

        if segment_gamma.size == 0:
            continue

        events.append(
            SeizureEvent(
                start_sec=start_sec,
                end_sec=end_sec,
                duration_sec=duration_sec,
                gamma_peak=float(np.max(segment_gamma)),
                pis_found=True,
                pis_start_sec=None if pis_s is None else float(t_full[pis_s]),
                pis_end_sec=None if pis_e is None else float(t_full[pis_e]),
            )
        )

    return events


# Backward-compatible alias for the old notebook function name.
detect_seizures_gamma_pis_demo = detect_seizures_gamma_pis
