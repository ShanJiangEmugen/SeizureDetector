import argparse
from pathlib import Path
import sys

# Allow running script directly from project root.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from seizure_detection.config import ensure_output_dir, load_config
from seizure_detection.data_io import scan_edf_files
from seizure_detection.batch import build_baseline_map, run_batch_detection
from seizure_detection.summary import build_animal_summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run EEG seizure detection from EDF files."
    )

    parser.add_argument(
        "--config",
        type=str,
        default="configs/default.yaml",
        help="Path to YAML config file.",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="Override input_dir in config.",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Override output_dir in config.",
    )

    parser.add_argument(
        "--channel",
        type=str,
        default=None,
        help="Override EEG channel in config.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    input_dir = args.input_dir or cfg.get("input_dir", "EDF")
    output_dir = args.output_dir or cfg.get("output_dir", "outputs")

    if args.channel is not None:
        channel = args.channel
    else:
        channel = cfg.get("channel", None)

    groups = cfg.get("groups", None)
    detector_kwargs = cfg.get("detector", {})

    baseline_cfg = cfg.get("baseline", {})
    baseline_method = baseline_cfg.get("method", "nth_edf")
    nth_index = int(baseline_cfg.get("nth_index", 1))

    output_dir = ensure_output_dir(output_dir)

    events_csv = output_dir / cfg.get("output_events_csv", "all_seizure_tracked.csv")
    summary_csv = output_dir / cfg.get("output_summary_csv", "animal_summary.csv")

    print("=" * 80)
    print("[INFO] Seizure detection started")
    print(f"[INFO] Config: {args.config}")
    print(f"[INFO] Input dir: {input_dir}")
    print(f"[INFO] Output dir: {output_dir}")
    print(f"[INFO] Channel: {channel}")
    print("=" * 80)

    edf_records = scan_edf_files(input_dir=input_dir, groups=groups)

    print(f"[INFO] Found {len(edf_records)} EDF files.")

    if len(edf_records) == 0:
        raise RuntimeError(f"No EDF files found under: {input_dir}")

    baseline_map = build_baseline_map(
        edf_records=edf_records,
        method=baseline_method,
        nth_index=nth_index,
    )

    print(f"[INFO] Built baseline map for {len(baseline_map)} animals.")

    df_events = run_batch_detection(
        edf_records=edf_records,
        baseline_map=baseline_map,
        channel=channel,
        detector_kwargs=detector_kwargs,
    )

    df_events.to_csv(events_csv, index=False)
    print(f"[DONE] Event table saved to: {events_csv}")

    df_summary = build_animal_summary(df_events)
    df_summary.to_csv(summary_csv, index=False)
    print(f"[DONE] Animal summary saved to: {summary_csv}")

    print("=" * 80)
    print("[DONE] Seizure detection finished")
    print("=" * 80)


if __name__ == "__main__":
    main()
