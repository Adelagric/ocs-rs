"""Join only complete, fingerprint-matched fresh operational measurements."""
import csv
import os
from pathlib import Path
from operational_data import finite, validate_algorithm, validate_native

base = Path(os.environ.get("REPRO_OUTPUT_DIR", "artifacts"))
rows = []
for panel in ("wheat", "mouse", "pig"):
    rpath = base / ("operational_ne_pig.csv" if panel == "pig" else "operational_ne_realpanels.csv")
    with rpath.open() as stream:
        measured = {float(r["target_ne"]): r for r in csv.DictReader(stream) if r["panel"] == panel}
    with (base / f"algorithm_operational_{panel}.csv").open() as stream:
        algorithm = {float(r["target_ne"]): r for r in csv.DictReader(stream)}
    if set(measured) != {25., 50., 100.} or set(algorithm) != set(measured):
        raise ValueError(f"incomplete operating-point measurements for {panel}")
    for target, r in measured.items():
        a = algorithm[target]
        validate_native(r)
        validate_algorithm(a)
        for key in ("ocsrs_s", "optisel_s"):
            finite(r, key, True)
        if float(r["ub"]) != float(a["ub"]) or any(r[k] != a[k] for k in ("K_md5", "bc_md5")):
            raise ValueError(f"different instances/caps in {panel}, target={target}")
        rows.append(dict(panel=panel, n=r["n"], m=r["m"], target_ne=target,
                         achieved_ne=r["achieved_ne"], proxy_delta_f_pct=r["proxy_delta_f_pct"],
                         support=r["support_ocsrs"], algorithm_s=a["algorithm_s"],
                         matrixfree_s=r["ocsrs_s"], optisel_s=r["optisel_s"],
                         gain_gap=r["gain_gap"], native_status=r["status"], algorithm_status=a["status"],
                         algorithm_verified=a["verified"], ub=r["ub"], K_md5=r["K_md5"], bc_md5=r["bc_md5"]))
with (base / "table2_operational.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
