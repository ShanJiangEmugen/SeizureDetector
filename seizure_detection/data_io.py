from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import mne
import numpy as np


def read_edf_channel(edf_path: str | Path, picks: Optional[str] = None):
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)

    if picks is None:
        picks = raw.ch_names[0]

    if picks not in raw.ch_names:
        raise ValueError(
            f"Channel '{picks}' not found in {edf_path}. "
            f"Available channels: {raw.ch_names}"
        )

    raw.pick([picks])

    x = raw.get_data()[0]
    fs = float(raw.info["sfreq"])

    return np.asarray(x, dtype=float), fs


def scan_edf_files(
    input_dir: str | Path,
    groups: Optional[Iterable[str]] = None,
) -> List[Tuple[str, str, Path]]:
    """
    Return list of (group, animal_id, edf_path).

    Expected structure:
    input_dir/
        A/
            Animal1/
                file.edf
        B/
            Animal2/
                file.edf
    """
    input_dir = Path(input_dir)

    if groups is None:
        groups = [p.name for p in input_dir.iterdir() if p.is_dir()]

    results = []

    for group in groups:
        group_dir = input_dir / group

        if not group_dir.exists():
            print(f"[WARN] Group folder not found: {group_dir}")
            continue

        for edf_path in sorted(group_dir.rglob("*.edf")):
            animal_id = edf_path.parent.name
            results.append((group, animal_id, edf_path))

    return results
