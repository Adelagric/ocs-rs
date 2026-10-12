"""Acceptance gates for fresh native/dense operating-point evidence."""
import math


def finite(row, key, positive=False):
    value = float(row[key])
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{key} must be finite" + (" and positive" if positive else ""))
    return value


def validate_native(row):
    if row["status"] != "Solved":
        raise ValueError("native operational result must have status Solved")
    target = finite(row, "target_ne", True)
    if target not in (25., 50., 100.):
        raise ValueError("unknown contribution-Ne target")
    achieved = finite(row, "achieved_ne", True)
    if abs(achieved - target) / target > .015:
        raise ValueError("native result did not reach the contribution-Ne target")
    finite(row, "ub", True)
    finite(row, "gain_ocsrs")
    return target, achieved


def validate_algorithm(row):
    if str(row["verified"]) != "True" or row["status"] != "Solved":
        raise ValueError("dense operational result lacks a verified Solved certificate")
    target = finite(row, "target_ne", True)
    achieved = finite(row, "achieved_ne", True)
    if target not in (25., 50., 100.) or abs(achieved - target) / target > .02:
        raise ValueError("dense result did not reach the contribution-Ne target")
    finite(row, "ub", True)
    finite(row, "algorithm_s", True)
    return target, achieved
