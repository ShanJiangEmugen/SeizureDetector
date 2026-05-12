from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from .data_io import read_edf_channel
from .detector import detect_seizures_gamma_pis


def build_baseline_map(
    edf_records: List[Tuple[str, str, Path]],
    method: str = "nth_edf",
    nth_index: int = 1,
) -> Dict[str, Path]:
    animal_to_paths = defaultdict(list)

    for group, animal_id, edf_path in edf_records:
        animal_to_paths[animal_id].append(edf_path)

    baseline_map = {}

    if method != "nth_edf":
        raise ValueError(f"Unsupported baseline method: {method}")

    for animal_id, paths in animal_to_paths.items():
        paths = sorted(paths)

        if len(paths) <= nth_index:
            print(
                f"[WARN] Animal {animal_id} has only {len(paths)} EDF file(s), "
                f"cannot use nth_index={nth_index} as baseline."
            )
            continue

        baseline_map[animal_id] = paths[nth_index]

    return baseline_map


def events_to_rows(events, group, animal_id, edf_path, baseline_path):
    rows = []

    for i, e in enumerate(events):
        rows.append(
            {
                "group": group,
                "animal_id": animal_id,
                "test_file": Path(edf_path).name,
                "test_path": str(edf_path),
                "baseline_file": Path(baseline_path).name,
                "baseline_path": str(baseline_path),
                "event_index": i,
                "start_sec": e.start_sec,
                "end_sec": e.end_sec,
                "duration_sec": e.duration_sec,
                "gamma_peak": e.gamma_peak,
                "pis_found": e.pis_found,
                "pis_start_sec": e.pis_start_sec,
                "pis_end_sec": e.pis_end_sec,
            }
        )

    return rows


def run_batch_detection(
    edf_records: List[Tuple[str, str, Path]],
    baseline_map: Dict[str, Path],
    channel: str | None = None,
    detector_kwargs: dict | None = None,
) -> pd.DataFrame:
    if detector_kwargs is None:
        detector_kwargs = {}

    all_rows = []

    for group, animal_id, edf_path in edf_records:
        baseline_path = baseline_map.get(animal_id)

        if baseline_path is None:
            print(f"[WARN] Skip {animal_id}: no baseline found.")
            continue

        print(f"[INFO] Processing animal={animal_id}, file={edf_path.name}")

        try:
            x, fs = read_edf_channel(edf_path, picks=channel)
            x_base, fs_base = read_edf_channel(baseline_path, picks=channel)

            if abs(fs - fs_base) > 1e-6:
                print(
                    f"[WARN] Sampling rate mismatch for {animal_id}: "
                    f"test fs={fs}, baseline fs={fs_base}"
                )

            events = detect_seizures_gamma_pis(
                x_full=x,
                fs=fs,
                x_baseline=x_base,
                **detector_kwargs,
            )

            rows = events_to_rows(
                events=events,
                group=group,
                animal_id=animal_id,
                edf_path=edf_path,
                baseline_path=baseline_path,
            )

            all_rows.extend(rows)

            print(f"[INFO] Found {len(events)} event(s).")

        except Exception as e:
            print(f"[ERROR] Failed on {edf_path}: {e}")

    return pd.DataFrame(all_rows)
