import os
import matplotlib.pyplot as plt
import numpy as np


def plot_seizure_timeline(t, signal, events, save_path=None):
    plt.figure(figsize=(12, 4))

    plt.plot(t, signal, linewidth=0.5)

    for e in events:
        plt.axvspan(e.start_sec, e.end_sec, alpha=0.3)

    plt.xlabel("Time (sec)")
    plt.ylabel("Signal")

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_each_seizure_segment(signal, fs, events, out_dir, prefix="seizure"):
    os.makedirs(out_dir, exist_ok=True)

    for i, e in enumerate(events):
        s = int(e.start_sec * fs)
        e_idx = int(e.end_sec * fs)

        seg = signal[s:e_idx]
        t = np.arange(len(seg)) / fs

        plt.figure(figsize=(6, 3))
        plt.plot(t, seg)
        plt.title(f"{prefix}_{i}")

        path = os.path.join(out_dir, f"{prefix}_{i}.png")
        plt.savefig(path, dpi=150)
        plt.close()
