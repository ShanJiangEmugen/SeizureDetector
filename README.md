# 🧠 Seizure-Detection-Toolkit

A modular EEG-based seizure detection pipeline with gamma-power analysis and post-ictal suppression validation.

<p align="left">
  <!-- Release -->
  <img src="https://img.shields.io/github/v/release/ShanJiangEmugen/seizure-detection-toolkit?color=blue&label=Release&style=flat-square" />
  <!-- License -->
  <img src="https://img.shields.io/github/license/ShanJiangEmugen/seizure-detection-toolkit?style=flat-square" />
  <!-- Issues -->
  <img src="https://img.shields.io/github/issues/ShanJiangEmugen/seizure-detection-toolkit?style=flat-square" />
  <!-- Python Version -->
  <img src="https://img.shields.io/badge/Python-3.8%20|%203.9%20|%203.10-blue?style=flat-square" />
  <!-- Platform -->
  <img src="https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey?style=flat-square" />
</p>

---

## Overview

**Seizure-Detection-Toolkit** is a lightweight, modular, and production-ready toolkit for **automated seizure detection from EEG recordings**.

The pipeline is designed for **rodent and primate electrophysiology experiments**, integrating:

- Gamma-band power analysis
- Adaptive artifact (super-peak) masking
- Statistical thresholding using baseline EEG
- Post-ictal suppression (PIS) validation
- Batch processing across multiple animals
- Structured outputs for downstream analysis

This toolkit is especially suitable for:

- Neuroscience research (epilepsy, ASD, gene therapy)
- Long-term EEG monitoring studies
- High-throughput behavioral / neural pipelines

---

## Features

### EEG-Based Seizure Detection
- Gamma-band power extraction via sliding window FFT
- Baseline-driven thresholding (mean + SD)
- Robust excursion detection (sustained high gamma)

### Artifact Handling
- Adaptive masking of large-amplitude peaks
- Prevents false positives from motion / noise artifacts
- Handles NaN interpolation in frequency domain

### PIS Validation (Key Feature)
- Detects post-ictal suppression-like patterns
- Uses dynamic ratio-based gamma drop
- Filters out non-seizure high-gamma events

### Batch Processing Pipeline
- Automatic EDF scanning (A/B/C groups)
- Animal-level baseline assignment
- Scalable across large datasets

### Clean Modular Codebase
- Separation of:
  - `detector` (core algorithm)
  - `data_io` (EDF handling)
  - `batch` (pipeline logic)
  - `summary` (aggregation)

---

## Installation

```bash
git clone https://github.com/your-username/seizure-detection-toolkit.git
pip install -r requirements.txt
```

## Quick Start

### 1. Prepare Data Structure

Expected directory format:

```text
EDF/
├── A/
│   ├── Animal1/
│   │   ├── file1.edf
│   │   ├── file2.edf
│   │   └── ...
├── B/
├── C/
cd seizure-detection-toolkit
```

Each animal should have multiple EDF files. By default, the **second EDF file** (sorted by filename) will be used as the baseline.

---

### 2. Run Seizure Detection

Basic command:
```
python scripts/run_batch_detection.py --config configs/default.yaml
```


This will:
- Scan all EDF files under the configured directory
- Assign baseline per animal
- Run seizure detection
- Save results to output folder

---

### 3. Optional: Override Parameters

You can override config settings directly from CLI:
```
python scripts/run_batch_detection.py
--config configs/default.yaml
--input-dir EDF
--output-dir outputs
--channel "EEG 2_FCXA-B"
```

Common options:

- `--input-dir`  
  Path to EDF dataset

- `--output-dir`  
  Folder where results will be saved

- `--channel`  
  Specific EEG channel name (default: first channel)

---

### 4. Output Files

After running, results will be saved to:

```text
outputs/
├── all_seizure_tracked.csv
└── animal_summary.csv
```


---

## Configuration (YAML)

Example configuration file: `configs/default.yaml`

```text
input_dir: "EDF"
output_dir: "outputs"

channel: null

groups:

"A"
"B"
"C"

baseline:
method: "nth_edf"
nth_index: 1

detector:
gamma_band: [20.0, 50.0]
win_sec: 1.0
step_sec: 1.0
thr_sd: 3.0

min_excursion_sec: 5.0

pis_search_sec: 90.0
min_pis_sec: 6.0
pis_ratio: 100.0

enable_peak_mask: true
peak_thr: 0.0015
```


---

## Output Format

### Event-Level Table (`all_seizure_tracked.csv`)

Each row represents one detected seizure event:

| column | description |
|------|------------|
| group | experimental group (A/B/C) |
| animal_id | subject ID |
| test_file | EDF filename |
| start_sec | seizure start time (seconds) |
| end_sec | seizure end time (seconds) |
| duration_sec | seizure duration |
| gamma_peak | peak gamma power during event |
| pis_start_sec | start of post-ictal suppression |
| pis_end_sec | end of post-ictal suppression |

---

### Animal Summary (`animal_summary.csv`)

Aggregated statistics per animal:

| column | description |
|------|------------|
| animal_id | subject ID |
| seizure_count | total number of detected seizures |
| total_duration_sec | total seizure duration |
| mean_duration_sec | average seizure duration |

---

## Example Workflow

1. Record EEG data
2. Export to EDF format
3. Organize files into A/B/C group folders
4. Run detection pipeline
5. Analyze:
   - seizure count
   - seizure duration
   - group-level differences

---

## Notes

- Designed for research use only
- Detection performance depends on:
  - EEG signal quality
  - sampling rate
  - parameter tuning
- Strongly recommended to validate results with visual inspection (QC plots will be added in future versions)






