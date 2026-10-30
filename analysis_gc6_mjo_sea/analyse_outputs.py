"""Assess existing diagnostic products without modifying their source files."""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from cf_units import Unit
from PIL import Image, ImageDraw
from scipy.interpolate import RegularGridInterpolator


HERE = Path(__file__).resolve().parent
ROOT = Path("/data/users/prince.xavier/SEAPy_climate")
RUNS = {
    "control": ("u-cv318", "u-cv318_GC5.0_CPL_C144"),
    "experiment": ("u-ea703", "u-ea703_GC6.0_CPL_C144"),
    "obs": ("obs", "obs_ERA5_GPM"),
}
GROUPS = ("MJOPy_output", "SEAPy_output")


def write_csv(name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def folder(group, role):
    return ROOT / group / RUNS[role][1]


def load_field(path):
    with netCDF4.Dataset(path) as dataset:
        candidates = [variable for variable in dataset.variables.values()
                      if variable.ndim >= 2 and not variable.name.endswith("bnds")
                      and "bnds" not in variable.dimensions]
        if len(candidates) != 1:
            raise ValueError(f"Ambiguous data variables: {path}")
        variable = candidates[0]
        values = np.ma.filled(variable[:].astype(float), np.nan)
        aliases = {"lat": "latitude", "lon": "longitude"}
        coords = {aliases.get(dim, dim): np.asarray(dataset.variables[dim][:], dtype=float)
                  for dim in variable.dimensions}
        return {"data": values, "dims": tuple(aliases.get(dim, dim) for dim in variable.dimensions), "coords": coords,
                "units": getattr(variable, "units", "unknown"), "name": variable.name}


def sample(field, latitudes, longitudes):
    if field["dims"][-2:] != ("latitude", "longitude"):
        raise ValueError(f"Unsupported dimensions: {field['dims']}")
    lat = field["coords"]["latitude"]
    lon = field["coords"]["longitude"]
    values = field["data"]
    lat_order = np.argsort(lat)
    lon_order = np.argsort(lon)
    values = values[..., lat_order, :][..., lon_order]
    interpolator = RegularGridInterpolator(
        (lat[lat_order], lon[lon_order]), np.moveaxis(values, (-2, -1), (0, 1)),
        bounds_error=False, fill_value=np.nan)
    grid_lat, grid_lon = np.meshgrid(latitudes, longitudes, indexing="ij")
    result = interpolator(np.column_stack((grid_lat.ravel(), grid_lon.ravel())))
    result = result.reshape((len(latitudes), len(longitudes)) + values.shape[:-2])
    return np.moveaxis(result, (0, 1), (-2, -1))


def weighted_mean(values, weights):
    valid = np.isfinite(values)
    return float(np.sum(np.where(valid, values, 0) * weights) /
                 np.sum(np.where(valid, weights, 0)))


def compare_arrays(arrays, weights):
    valid = np.logical_and.reduce([np.isfinite(array) for array in arrays])
    weights = np.broadcast_to(weights, arrays[0].shape)
    if not np.any(valid):
        raise ValueError("No common valid cells")
    arrays = [np.where(valid, array, np.nan) for array in arrays]
    result = {"n_common": int(valid.sum()), "obs_mean": weighted_mean(arrays[2], weights)}
    for role, values in zip(("control", "experiment"), arrays[:2]):
        difference = values - arrays[2]
        result[f"{role}_mean"] = weighted_mean(values, weights)
        result[f"{role}_bias"] = weighted_mean(difference, weights)
        result[f"{role}_rmse"] = np.sqrt(weighted_mean(difference ** 2, weights))
        centred = values - weighted_mean(values, weights)
        reference = arrays[2] - result["obs_mean"]
        denominator = np.sqrt(weighted_mean(centred ** 2, weights) *
                              weighted_mean(reference ** 2, weights))
        result[f"{role}_correlation"] = (
            weighted_mean(centred * reference, weights) / denominator if denominator else np.nan)
        result[f"{role}_rms"] = np.sqrt(weighted_mean(values ** 2, weights))
    result["rmse_change_percent"] = (
        100 * (result["experiment_rmse"] / result["control_rmse"] - 1)
        if result["control_rmse"] else np.nan)
    result["verdict"] = "improved" if result["experiment_rmse"] < result["control_rmse"] else "deteriorated"
    return result


def scalar_metrics():
    rows = []
    conflicts = []
    for group in GROUPS:
        tables = {}
        for role, (runid, _) in RUNS.items():
            primary = folder(group, role) / f"{runid}_metrics.csv"
            with primary.open() as handle:
                tables[role] = {key: float(value) for key, value in csv.reader(handle)}
            alternative = folder(group, role) / "metrics.csv"
            if alternative.exists():
                with alternative.open() as handle:
                    secondary = {key: float(value) for key, value in csv.reader(handle)}
                for key in tables[role].keys() & secondary.keys():
                    if tables[role][key] != secondary[key]:
                        conflicts.append({"group": group, "role": role, "metric": key,
                                          "primary": tables[role][key], "generic": secondary[key]})
        for key in sorted(set.intersection(*(set(table) for table in tables.values()))):
            row = {"group": group, "metric": key,
                   **{role: tables[role][key] for role in RUNS}}
            control_error = abs(row["control"] - row["obs"])
            experiment_error = abs(row["experiment"] - row["obs"])
            row.update(control_absolute_error=control_error, experiment_absolute_error=experiment_error,
                       error_change=experiment_error - control_error,
                       verdict="unchanged" if np.isclose(control_error, experiment_error) else
                       "improved" if experiment_error < control_error else "deteriorated")
            rows.append(row)
    write_csv("scalar_metrics.csv", rows)
    return conflicts


def field_metrics():
    rows = []
    skipped = []
    for group in GROUPS:
        for path in sorted(folder(group, "obs").glob("obs_*.nc")):
            suffix = path.name[4:]
            if suffix.endswith("_filtered.nc"):
                continue
            paths = [folder(group, role) / f"{runid}_{suffix}" for role, (runid, _) in RUNS.items()]
            if not all(item.exists() for item in paths):
                skipped.append({"group": group, "product": suffix, "reason": "missing matched product"})
                continue
            fields = [load_field(item) for item in paths]
            if len({field["dims"] for field in fields}) != 1:
                skipped.append({"group": group, "product": suffix,
                                "reason": f"dimension mismatch: {[field['dims'] for field in fields]}"})
                continue
            if len({field["data"].shape[:-2] for field in fields}) != 1:
                skipped.append({"group": group, "product": suffix,
                                "reason": "nonspatial shape mismatch"})
                continue
            units = [field["units"] for field in fields]
            if len(set(units)) > 1:
                allowed = {"kg m-2 day-1", "mm day-1", "mm/day", "kg m-2 d-1"}
                if not set(units).issubset(allowed):
                    target = Unit(units[0])
                    if all(Unit(unit).is_convertible(target) for unit in units):
                        for field in fields:
                            field["data"] = Unit(field["units"]).convert(field["data"], target)
                    else:
                        skipped.append({"group": group, "product": suffix,
                                        "reason": f"incompatible unit metadata: {units}"})
                        continue
            if fields[0]["dims"] == ("lead_lag", "longitude"):
                assert all(np.array_equal(fields[0]["coords"][dim], field["coords"][dim])
                           for dim in fields[0]["dims"] for field in fields[1:])
                arrays = [field["data"] for field in fields]
                row = {"group": group, "product": suffix, "domain": "40-180E, all saved lags", "units": units[0]}
                rows.append({**row, **compare_arrays(arrays, np.ones_like(arrays[0]))})
                continue
            if fields[0]["dims"][-2:] != ("latitude", "longitude"):
                skipped.append({"group": group, "product": suffix, "reason": "spectral coordinates; scalar metrics used"})
                continue
            if group == "SEAPy_output" and "variance" in suffix:
                domains = [("SEA wave box", 100, 130, -10, 10)]
            elif "percentile" in suffix or group == "SEAPy_output":
                domains = [("SEA", 95, 130, -10, 25)]
            else:
                domains = [("Indo-Pacific", 40, 180, -15, 15), ("SEA", 95, 130, -10, 25)]
            for domain, west, east, south, north in domains:
                latitudes = np.arange(south, north + 0.1, 2.5)
                longitudes = np.arange(west, east + 0.1, 2.5)
                arrays = [sample(field, latitudes, longitudes) for field in fields]
                weights = np.cos(np.deg2rad(latitudes))[:, None]
                row = {"group": group, "product": suffix, "domain": domain, "units": units[0]}
                rows.append({**row, **compare_arrays(arrays, weights)})
            if group == "SEAPy_output" and "variance" in suffix:
                row = {"group": group, "product": suffix, "domain": "native arithmetic 100-130E,10S-12N"}
                for role, field in zip(RUNS, fields):
                    lat_mask = (field["coords"]["latitude"] >= -10) & (field["coords"]["latitude"] <= 12)
                    lon_mask = (field["coords"]["longitude"] >= 100) & (field["coords"]["longitude"] <= 130)
                    row[f"{role}_mean"] = float(np.nanmean(field["data"][..., lat_mask, :][..., lon_mask]))
                rows.append(row)
    write_csv("field_metrics.csv", rows)
    write_csv("excluded_products.csv", skipped)
    return rows


def dates_and_indices():
    result = {}
    for role, (runid, _) in RUNS.items():
        rmm = np.loadtxt(folder("MJOPy_output", role) / f"RMMs_{runid}.txt")
        valid = (rmm[:, 6] >= 0) & (rmm[:, 5] >= 1) & (rmm[:, 5] <= 8)
        active = valid & (rmm[:, 6] >= 1)
        stats = {"rmm_first_date": rmm[0, :3].astype(int).tolist(),
                 "rmm_last_date": rmm[-1, :3].astype(int).tolist(),
                 "rmm_days": len(rmm), "rmm_valid_days": int(valid.sum()),
                 "rmm_active_percent": float(100 * active.sum() / valid.sum()),
                 "rmm_active_phase_counts": np.bincount(rmm[active, 5].astype(int), minlength=9)[1:].tolist(),
                 "rmm_calendar_valid_first": rmm[valid][0, :3].astype(int).tolist(),
                 "rmm_calendar_valid_last": rmm[valid][-1, :3].astype(int).tolist()}
        for group in GROUPS:
            stats[group] = {}
            for path in sorted(folder(group, role).glob(f"{runid}_*dates.txt")):
                values = np.loadtxt(path, dtype=int, ndmin=1)
                months = values // 100 % 100
                stats[group][path.name] = {"count": len(values), "first": int(values.min()) if len(values) else None,
                                          "last": int(values.max()) if len(values) else None,
                                          "outside_NDJF": int(((months > 2) & (months < 11)).sum())}
        result[role] = stats
    (HERE / "sampling_and_events.json").write_text(json.dumps(result, indent=2))
    return result


def inventory():
    rows = []
    for group in GROUPS:
        for role in RUNS:
            for path in sorted(folder(group, role).glob("*")):
                if path.is_file():
                    stat = path.stat()
                    rows.append({"group": group, "role": role, "path": str(path),
                                 "bytes": stat.st_size, "mtime_epoch": stat.st_mtime})
    write_csv("source_inventory.csv", rows)


def contact_sheet(name, group, suffixes):
    width, height = 750, 520
    canvas = Image.new("RGB", (width * 3, height * len(suffixes)), "white")
    draw = ImageDraw.Draw(canvas)
    for row_index, suffix in enumerate(suffixes):
        for column, (role, (runid, _)) in enumerate(RUNS.items()):
            path = folder(group, role) / f"{runid}_{suffix}.png"
            if not path.exists():
                continue
            with Image.open(path) as image:
                image = image.convert("RGB")
                image.thumbnail((width - 10, height - 32))
                canvas.paste(image, (column * width + (width - image.width) // 2, row_index * height + 27))
            draw.text((column * width + 10, row_index * height + 8), f"{role}: {runid} | {suffix}", fill="black")
    canvas.save(HERE / name)


def figures(rows):
    plt.rcParams.update({"font.size": 9, "savefig.dpi": 180})
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.3), constrained_layout=True)
    suffixes = ["precipitation_flux_LeadLagCorr.nc", "toa_outgoing_longwave_flux_LeadLagCorr.nc"]
    for row_index, suffix in enumerate(suffixes):
        for column, (role, (runid, _)) in enumerate(RUNS.items()):
            field = load_field(folder("MJOPy_output", role) / f"{runid}_{suffix}")
            ax = axes[row_index, column]
            artist = ax.contourf(field["coords"]["longitude"], field["coords"]["lead_lag"],
                                 field["data"], levels=np.linspace(-0.9, 0.9, 19), cmap="RdBu_r", extend="both")
            ax.set_title(f"{role.title()} ({runid})")
            ax.set_xlabel("Longitude (degrees E)")
            ax.set_ylabel("Lag (days)")
            ax.grid(alpha=0.15)
    fig.colorbar(artist, ax=axes, shrink=0.75, label="Saved lag-correlation coefficient")
    fig.savefig(HERE / "figure1_mjo_propagation.png")
    plt.close(fig)
    contact_sheet("inspection_mjo.png", "MJOPy_output", [
        "precipitation_flux_Mean_winter", "precipitation_flux_Filt_Var_winter",
        "precipitation_flux_wavenum_freq_season_winter", "precipitation_flux_winter_phase_composites",
        "phase_freq_plot", "precipitation_flux_filt_prop_comp"])
    contact_sheet("inspection_sea.png", "SEAPy_output", [
        "CS_precip_winds850_composite", "CES_precip_winds850_composite",
        "CS_sst_winds850_composite", "precipitation_flux_MJO_ndjf_stddev",
        "precipitation_flux_Kelvin_jjas_stddev", "NDJF_precip_95_percentile",
        "CS_precip_95_percentile_fraction"])
    print("Field comparisons:", len(rows))


def self_test():
    control = np.array([[2., 2.], [2., np.nan]])
    experiment = np.array([[1., 1.], [1., 100.]])
    obs = np.zeros((2, 2))
    result = compare_arrays([control, experiment, obs], np.ones((2, 1)))
    assert result["control_rmse"] == 2
    assert result["experiment_rmse"] == 1
    assert result["n_common"] == 3
    assert result["rmse_change_percent"] == -50
    field = {"dims": ("latitude", "longitude"), "coords": {
        "latitude": np.array([0., 2.]), "longitude": np.array([0., 2.])},
        "data": np.array([[0., 2.], [2., 4.]])}
    assert np.allclose(sample(field, np.array([1.]), np.array([1.])), [[2.]])
    print("PASS: common masking, RMSE direction and bilinear interpolation")


if __name__ == "__main__":
    import sys
    self_test()
    if "--self-test" not in sys.argv:
        inventory()
        conflicts = scalar_metrics()
        (HERE / "metric_conflicts.json").write_text(json.dumps(conflicts, indent=2))
        print("Metric conflicts:", conflicts)
        print(json.dumps(dates_and_indices(), indent=2))
        figures(field_metrics())