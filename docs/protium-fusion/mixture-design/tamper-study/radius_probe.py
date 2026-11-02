#!/usr/bin/env python3
"""Rerun an untamped, prescribed-inflow N14 pulse at a chosen radius."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "b2-completion"))
from radial_completion import simulate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("radius_m", type=float)
    parser.add_argument("--output", type=Path, required=True,
                        help="Output prefix; .json and .npz will be written")
    args = parser.parse_args()
    if args.radius_m <= 0:
        parser.error("radius_m must be positive")
    profile = ROOT / "inputs/pulse-profiles/scaled-hot-A-64-R30-initial.json"
    initial = json.loads(profile.read_text())
    initial["outer_radius_m"] = (
        np.asarray(initial["outer_radius_m"]) * args.radius_m / 30
    ).tolist()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = simulate(
        "A", 3e8, 1e8, args.radius_m, zones=64, nuclear=True,
        reduced_a=False, radiation=True, kappa=0.001, initial_profile=initial,
        final_time_units=4, rtol=1e-6, atol=1e-11,
        viscosity_mode="compression_cap", trace_samples=60,
        save_trace=args.output.with_suffix(".npz"),
    )
    args.output.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result.get(key) for key in (
        "success", "termination_reason", "mass_kg", "initial_energy_J",
        "nuclear_energy_J", "deposited_nuclear_over_initial_energy",
        "desired_product_per_initial_catalyst", "peak_temperature_K",
        "radiation_escaped_J", "energy_residual_fraction_initial",
    )}, indent=2))


if __name__ == "__main__":
    main()
