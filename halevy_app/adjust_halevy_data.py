"""Reproduce rounding-adjusted Halevy data.

This script orthogonally projects each observed bundle (X, Y) onto its
budget line defined by the reported intercepts, snapping to the nearest
endpoint if the projection lies outside the feasible segment. The
adjustment is deterministic and minimizes Euclidean distance while
ensuring the implied income is exactly 1.

Defaults use the repository CSV locations, but input/output paths can be
overridden via command-line arguments.
"""

import argparse
import os
from pathlib import Path

# Set conservative threading to avoid SHM errors under sandboxing.
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("KMP_INIT_AT_FORK", "FALSE")

import numpy as np
import pandas as pd


def adjust(df: pd.DataFrame) -> pd.DataFrame:
    required = ["Subject", "Observation", "X", "Y", "X-intercept", "Y-intercept"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    px = 1.0 / df["X-intercept"].to_numpy()
    py = 1.0 / df["Y-intercept"].to_numpy()
    X = df["X"].to_numpy()
    Y = df["Y"].to_numpy()

    income_before = px * X + py * Y

    # Orthogonal projection of (X, Y) onto p·q = 1.
    den = px**2 + py**2
    lambda_factor = (income_before - 1.0) / den
    X_proj = X - lambda_factor * px
    Y_proj = Y - lambda_factor * py

    x_int = df["X-intercept"].to_numpy()
    y_int = df["Y-intercept"].to_numpy()

    inside = (X_proj >= 0) & (X_proj <= x_int) & (Y_proj >= 0) & (Y_proj <= y_int)

    X_adj = X_proj.copy()
    Y_adj = Y_proj.copy()

    need_clip = ~inside
    if need_clip.any():
        # Distances to endpoints (0, y_int) and (x_int, 0)
        d1 = (X - 0.0) ** 2 + (Y - y_int) ** 2
        d2 = (X - x_int) ** 2 + (Y - 0.0) ** 2

        use_first = d1 <= d2
        idx1 = need_clip & use_first
        idx2 = need_clip & ~use_first

        X_adj[idx1] = 0.0
        Y_adj[idx1] = y_int[idx1]

        X_adj[idx2] = x_int[idx2]
        Y_adj[idx2] = 0.0

    income_after = px * X_adj + py * Y_adj

    out_df = df.copy()
    out_df["X_old"] = out_df["X"]
    out_df["Y_old"] = out_df["Y"]
    out_df["X"] = X_adj
    out_df["Y"] = Y_adj
    out_df["income_before"] = income_before
    out_df["income_after"] = income_after
    out_df["abs_delta"] = np.hypot(out_df["X"] - out_df["X_old"], out_df["Y"] - out_df["Y_old"])

    return out_df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data") / "Halevy et al (2016) - Data.csv",
        help="Path to source CSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data") / "Halevy_Data_adjusted.csv",
        help="Where to write adjusted CSV",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    out_df = adjust(df)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, index=False)

    print(f"Wrote {len(out_df)} rows to {args.output}")
    print(f"max abs_delta={out_df['abs_delta'].max():.10f}")
    print(f"mean abs_delta={out_df['abs_delta'].mean():.10f}")
    print(
        f"income_after range=[{out_df['income_after'].min():.16f}, {out_df['income_after'].max():.16f}]"
    )


if __name__ == "__main__":
    main()
