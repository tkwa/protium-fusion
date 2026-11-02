#!/usr/bin/env python3
"""Reproduce table interpolation and ideal assembly bounds; standard library only."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def interpolate(rows, x):
    assert rows[0][0] <= x <= rows[-1][0]
    for (x0, y0), (x1, y1) in zip(rows, rows[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    raise AssertionError("uncovered interpolation")


def calculate():
    data = json.loads((HERE / "data.json").read_text())
    tab = data["roth_compression_annealed_basal"]
    inp = data["assembly_screen_inputs"]
    density = tab["density_kg_m3"]
    fill = inp["archived_bulk_density_kg_m3"] / density
    points = []
    for temp in [1000, 1100, 1200, 1250]:
        k = interpolate(tab["rows"], temp - 273.15)
        points.append({"T_K": temp, "material_basal_k_W_m_K": k,
                       "parallel_packing_upper_k_W_m_K": fill * k,
                       "material_to_required1250_ratio": k / inp["plant_required_constant_k_for_1250K"],
                       "parallel_packing_to_required1250_ratio": fill * k / inp["plant_required_constant_k_for_1250K"]})
    sigma = inp["mass_per_area_kg_m2"]
    q = inp["peak_heat_flux_W_m2"]
    result = {"schema": "graphite-material-geometry-screen-v1",
              "input_sha256": hashlib.sha256((HERE / "data.json").read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "interpolation": "Piecewise linear in source temperature; no extrapolation or uncertainty distribution.",
              "points": points,
              "density_and_geometry": {
                  "solid_fill_fraction_for_archived_density": fill,
                  "dense_material_thickness_m": sigma / density,
                  "archived_density_thickness_m": sigma / inp["archived_bulk_density_kg_m3"],
                  "packing_bound_scope": "Voigt parallel paths with zero void conductivity and perfect contacts. It is an upper bound for this assumed packing, not a measured assembly."},
              "contact_budget": [{"allocated_drop_K": d, "maximum_added_R_m2_K_W": d / q}
                                 for d in [5, 10, 20]],
              "existing_surface_R_m2_K_W": inp["surface_laminate_thickness_m"] / inp["surface_laminate_k_W_m_K"] + inp["post_resistance_m2_K_W"],
              "required1200_to_material1200_ratio": inp["plant_required_constant_k_for_1200K"] / points[2]["material_basal_k_W_m_K"]}
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    check("Celsius Kelvin conversion", abs((750.9 + 273.15) - 1024.05) < 1e-10)
    check("interpolation preserves exact table knot", interpolate(tab["rows"], 750.9) == 404.206)
    check("positive and ordered source temperatures", all(x[0] < y[0] and x[1] > 0 for x,y in zip(tab["rows"],tab["rows"][1:])))
    check("fill fraction bounded", 0 < fill <= 1)
    for row in points:
        check(f"packing no higher than material at{row['T_K']}K", 0 < row["parallel_packing_upper_k_W_m_K"] <= row["material_basal_k_W_m_K"])
    check("constant325 material attainable over reviewed window", min(p["material_basal_k_W_m_K"] for p in points) > inp["plant_required_constant_k_for_1250K"])
    check("constant934 above measured material at1200K", inp["plant_required_constant_k_for_1200K"] > points[2]["material_basal_k_W_m_K"])
    check("dense thickness mass identity", abs(density * result["density_and_geometry"]["dense_material_thickness_m"] - sigma) < 1e-10)
    check("surface resistance unchanged", abs(result["existing_surface_R_m2_K_W"] - 4.15e-5) < 1e-15)
    for row in result["contact_budget"]:
        check(f"contact drop{row['allocated_drop_K']}K", abs(row["maximum_added_R_m2_K_W"] * q - row["allocated_drop_K"]) < 1e-12)
    result["checks"] = {"passed": len(checks), "names": checks}
    return result


if __name__ == "__main__":
    result = calculate()
    (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"checks": result["checks"]["passed"], "points": result["points"]}, indent=2))
