"""Small regression checks for reproduction integrity; no benchmark timings."""
import csv
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import runpy
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


class ReproductionTests(unittest.TestCase):
    def test_failed_figure_restore_retains_recovery_copy(self):
        script = (ROOT / "research/repro/repro.sh").read_text()
        function = "promote_outputs() {" + script.split("promote_outputs() {", 1)[1].split("\nrun_stage()", 1)[0]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for subdir in ("research", "fresh", "tools"):
                (root / subdir).mkdir()
            for suffix in ("pdf", "png"):
                (root / f"research/fig_frontier.{suffix}").write_text(f"reported {suffix}")
                (root / f"fresh/fig_frontier.{suffix}").write_text(f"fresh {suffix}")
            mv = root / "tools/mv"
            mv.write_text('#!/bin/sh\ncase "$1" in *.previous.tmp|*.png.repro-*.tmp) exit 7;; esac\nexec /bin/mv "$@"\n')
            mv.chmod(0o755)
            env = dict(os.environ, ROOT=str(root), REPRO_FIGURE_DIR=str(root / "fresh"),
                       REPRO_UPDATE_FIGURES="1", PATH=str(root / "tools") + ":/usr/bin:/bin")
            result = subprocess.run(["/bin/bash", "-c", function + '\npromote_outputs "fig_frontier.pdf fig_frontier.png"'],
                                    env=env, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            backups = list((root / "research").glob("fig_frontier.pdf.*.previous.tmp"))
            self.assertEqual(len(backups), 1, "erased the original figure after failed restoration")
            self.assertEqual(backups[0].read_text(), "reported pdf")

    def test_figures_reject_nonfinite_population_or_timing_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "research").mkdir()
            (root / "artifacts").mkdir()
            for name in ("fig_frontier.py", "fig_multigen.py"):
                shutil.copy(ROOT / "research" / name, root / "research" / name)
            with (root / "artifacts/table2_operational.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["panel", "target_ne", "proxy_delta_f_pct", "support", "matrixfree_s", "native_status", "algorithm_status", "algorithm_verified"])
                writer.writeheader()
                for panel in ("wheat", "mouse", "pig"):
                    for target in (25, 50, 100):
                        writer.writerow(dict(panel=panel, target_ne=target, proxy_delta_f_pct=50 / target,
                                             support="nan", matrixfree_s="nan", native_status="Solved",
                                             algorithm_status="Solved", algorithm_verified="True"))
            source = ROOT / "research/repro/reported/v0.4.0/multigen_validation.csv"
            with source.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["gain"] = "nan"
            with (root / "artifacts/multigen_validation.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            for name in ("fig_frontier.py", "fig_multigen.py"):
                result = subprocess.run([sys.executable, "research/" + name], cwd=root,
                                        capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0, name + " accepted NaN")
                self.assertIn("finite", result.stderr)

    def test_failed_second_figure_copy_preserves_both_reported_formats(self):
        script = (ROOT / "research/repro/repro.sh").read_text()
        function = "promote_outputs() {" + script.split("promote_outputs() {", 1)[1].split("\nrun_stage()", 1)[0]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "research").mkdir()
            (root / "fresh").mkdir()
            (root / "tools").mkdir()
            for suffix in ("pdf", "png"):
                (root / f"research/fig_frontier.{suffix}").write_text(f"reported {suffix}")
                (root / f"fresh/fig_frontier.{suffix}").write_text(f"fresh {suffix}")
            cp = root / "tools/cp"
            cp.write_text('#!/bin/sh\ncase "$1" in *.png) exit 7;; esac\nexec /bin/cp "$@"\n')
            cp.chmod(0o755)
            env = dict(os.environ, ROOT=str(root), REPRO_FIGURE_DIR=str(root / "fresh"),
                       REPRO_UPDATE_FIGURES="1", PATH=str(root / "tools") + ":/usr/bin:/bin")
            result = subprocess.run(["/bin/bash", "-c", function + '\npromote_outputs "fig_frontier.pdf fig_frontier.png"'],
                                    env=env, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            for suffix in ("pdf", "png"):
                self.assertEqual((root / f"research/fig_frontier.{suffix}").read_text(), f"reported {suffix}")

    def test_failed_output_copy_is_a_stage_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(ROOT / "research", root / "research")
            tools = root / "tools"
            tools.mkdir()
            for name in ("dirname", "mkdir", "date", "rm", "cat"):
                (tools / name).symlink_to(shutil.which(name))
            for name, script in (("cargo", '#!/bin/sh\nprintf "fixture\\n" > "$REPRO_OUTPUT_DIR/scaling_matrixfree.csv"\n'),
                                 ("cp", "#!/bin/sh\nexit 7\n")):
                path = tools / name
                path.write_text(script)
                path.chmod(0o755)
            env = dict(os.environ, PATH=str(tools), REPRO_STEPS="scaling_matrixfree")
            result = subprocess.run(["/bin/bash", str(root / "research/repro/repro.sh")],
                                    env=env, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("FAILED scaling_matrixfree", result.stdout)

    def test_dense_cap_must_be_finite(self):
        import numpy as np
        from research.repro.dense_support import ConvergenceError, verify_solution
        for cap in (float("nan"), float("inf")):
            with self.assertRaisesRegex(ConvergenceError, "finite"):
                verify_solution(np.eye(4), np.zeros(4), np.array([True, False, True, False]),
                                np.full(4, .25), cap, np.zeros(2), 0.)

    def test_operational_cli_requires_panels_and_valid_native_evidence(self):
        result = subprocess.run([sys.executable, str(ROOT / "research/repro/sf_operational_ne.py")],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0, "no-panel invocation silently succeeded")
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            inputs = base / "inputs"
            inputs.mkdir()
            kp, bp = inputs / "wheat_K.csv", inputs / "wheat_bc.csv"
            kp.write_text("1,0,0,0\n0,1,0,0\n0,0,1,0\n0,0,0,1\n")
            bp.write_text("bv,oc,sex\n4,0,male\n3,0,female\n2,0,male\n1,0,female\n")
            with (base / "operational_ne_realpanels.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["panel", "target_ne", "achieved_ne", "ub", "gain_ocsrs", "status", "K_md5", "bc_md5"])
                writer.writeheader()
                for target in (25, 50, 100):
                    writer.writerow(dict(panel="wheat", target_ne=target, achieved_ne="nan", ub=.5,
                                         gain_ocsrs="nan", status="IterationLimit",
                                         K_md5=hashlib.md5(kp.read_bytes()).hexdigest(),
                                         bc_md5=hashlib.md5(bp.read_bytes()).hexdigest()))
            env = dict(os.environ, REPRO_OUTPUT_DIR=str(base), REPRO_INPUT_DIR=str(inputs))
            result = subprocess.run([sys.executable, str(ROOT / "research/repro/sf_operational_ne.py"), "wheat"],
                                    env=env, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertNotIn("verified=True", result.stdout)

    def test_scaling_rejects_infinite_build_time(self):
        loader = runpy.run_path(str(ROOT / "research/fig_scaling.py"))["load_data"]
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "scaling_matrixfree.csv").write_text(
                "n,m,g_build_s,mfree_solve_s,support,status\n1000,1000,inf,.1,10,Solved\n")
            (base / "scaling_structured.csv").write_text(
                "n,m,target_ne,ne,support,solve_s,status\n1000,1000,50,50,90,.2,Solved\n")
            with self.assertRaisesRegex(ValueError, "finite"):
                loader(base)

    def test_dense_input_rejects_indefinite_kinship(self):
        from research.repro.dense_support import load_panel
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "bench_K_4.csv").write_text("-1,0,0,0\n0,1,0,0\n0,0,1,0\n0,0,0,1\n")
            (base / "bench_bc_4.csv").write_text("bv,oc,sex\n4,0,male\n3,0,female\n2,0,male\n1,0,female\n")
            with patch.dict(os.environ, {"REPRO_LEGACY_DIR": str(base)}):
                with self.assertRaisesRegex(ValueError, "positive definite"):
                    load_panel("4", legacy=True)

    def test_scaling_rejects_nonfinite_measurements(self):
        loader = runpy.run_path(str(ROOT / "research/fig_scaling.py"))["load_data"]
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "scaling_matrixfree.csv").write_text(
                "n,m,g_build_s,mfree_solve_s,support,status\n1000,1000,0.1,0.1,10,Solved\n")
            template = "n,m,target_ne,ne,support,solve_s,status\n1000,1000,50,{ne},90,{time},Solved\n"
            for ne, time in (("nan", ".2"), ("50", "nan")):
                (base / "scaling_structured.csv").write_text(template.format(ne=ne, time=time))
                with self.assertRaisesRegex(ValueError, "finite"):
                    loader(base)

    def test_invalid_stage_names_and_empty_tokens_fail_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(ROOT / "research", root / "research")
            tools = root / "tools"
            tools.mkdir()
            for name in ("dirname", "mkdir", "date"):
                (tools / name).symlink_to(shutil.which(name))
            for selection in ("unknown", "table1,,structured", "table1,"):
                env = dict(os.environ, PATH=str(tools), REPRO_STEPS=selection)
                result = subprocess.run(["/bin/bash", str(root / "research/repro/repro.sh")],
                                        env=env, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0, selection)
                self.assertIn("Unknown stage:", result.stderr)

    def test_strict_missing_dependency_returns_omission_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(ROOT / "research", root / "research")
            tools = root / "tools"
            tools.mkdir()
            for name in ("dirname", "mkdir", "date", "rm"):
                (tools / name).symlink_to(shutil.which(name))
            env = dict(os.environ, PATH=str(tools), REPRO_STEPS="structured", REPRO_STRICT="1")
            result = subprocess.run(["/bin/bash", str(root / "research/repro/repro.sh")],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stdout)
            self.assertIn("Actual failures: 0; dependency/data omissions: 1", result.stdout)

    def test_dense_verified_solution_and_iteration_exhaustion(self):
        import numpy as np
        from research.repro.dense_support import ConvergenceError, support_first, verify_solution
        K = np.eye(4)
        bv = np.array([4., 3., 2., 1.])
        male = np.array([True, False, True, False])
        c, support, iterations = support_first(K, bv, male, .3)
        offset = np.sqrt(.05 / 4)
        np.testing.assert_allclose(c, [.25 + offset, .25 + offset, .25 - offset, .25 - offset])
        self.assertEqual(len(support), 4)
        self.assertGreater(iterations, 0)
        with self.assertRaisesRegex(ConvergenceError, "iteration limit"):
            support_first(K, bv, male, .3, max_iter=0)
        with self.assertRaisesRegex(ConvergenceError, "KKT"):
            verify_solution(K, bv, male, np.full(4, .25), .3, np.array([0., 0.]), 0.)

    def test_scaling_rejects_mismatch_and_accepts_no_projection_rows(self):
        loader = runpy.run_path(str(ROOT / "research/fig_scaling.py"))["load_data"]
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "scaling_matrixfree.csv").write_text(
                "n,m,g_build_s,mfree_solve_s,support,status\n1000,1000,0.1,0.1,10,Solved\n")
            source = base / "scaling_structured.csv"
            template = "n,m,target_ne,ne,support,solve_s,status\n1000,{m},{target},50,90,.2,Solved\n"
            source.write_text(template.format(m=1000, target=50))
            built, omitted, _, m = loader(base)
            self.assertEqual((len(built), len(omitted), m), (1, 0, 1000))
            source.write_text(template.format(m=500, target=50))
            with self.assertRaisesRegex(ValueError, "marker count"):
                loader(base)
            source.write_text(template.format(m=1000, target=100))
            with self.assertRaisesRegex(ValueError, "Ne_c=50"):
                loader(base)

    def test_omitted_figure_preserves_reported_manuscript_image(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(ROOT / "research", root / "research")
            tools = root / "tools"
            tools.mkdir()
            for name in ("dirname", "mkdir", "date", "rm"):
                (tools / name).symlink_to(shutil.which(name))
            figure = root / "research/fig_frontier.pdf"
            figure.write_bytes(b"reported figure")
            env = dict(os.environ, PATH=str(tools), REPRO_STEPS="fig_frontier")
            result = subprocess.run(["/bin/bash", str(root / "research/repro/repro.sh")],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertTrue(figure.exists(), "omission deleted reported manuscript image")
            self.assertEqual(figure.read_bytes(), b"reported figure")

    def test_failed_cargo_stage_returns_failure_and_removes_stale_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(ROOT / "research", root / "research")
            tools = root / "tools"
            tools.mkdir()
            for name in ("dirname", "mkdir", "date", "cp", "rm", "cat", "tee"):
                target = shutil.which(name)
                if target:
                    (tools / name).symlink_to(target)
            cargo = tools / "cargo"
            cargo.write_text("#!/bin/sh\nexit 7\n")
            cargo.chmod(0o755)
            (root / "artifacts").mkdir()
            stale = root / "artifacts" / "scaling_matrixfree.csv"
            stale.write_text("stale\n")
            env = dict(os.environ, PATH=str(tools), REPRO_STEPS="scaling_matrixfree")
            result = subprocess.run(["/bin/bash", str(root / "research/repro/repro.sh")],
                                    env=env, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertFalse(stale.exists(), "failed stage retained its old CSV")

    def test_scaling_refuses_missing_structured_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "research").mkdir()
            (root / "artifacts").mkdir()
            shutil.copy(ROOT / "research/fig_scaling.py", root / "research/fig_scaling.py")
            (root / "artifacts/scaling_matrixfree.csv").write_text(
                "n,m,g_build_s,mfree_solve_s,support,status\n"
                "1000,1000,0.1,0.1,10,Solved\n2000,1000,infeasible,0.2,10,Solved\n")
            result = subprocess.run([sys.executable, "research/fig_scaling.py"], cwd=root,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, "silently substituted loose-cap data")
            self.assertIn("scaling_structured.csv", result.stderr)

    def test_dense_prototype_rejects_infeasible_cap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shutil.copytree(ROOT / "research/repro", root / "repro")
            script = root / "repro/sf_at_ub.py"
            script.write_text(script.read_text().replace("/tmp/bench_", str(root / "bench_")))
            (root / "bench_K_4.csv").write_text("1,0,0,0\n0,1,0,0\n0,0,1,0\n0,0,0,1\n")
            (root / "bench_bc_4.csv").write_text("bv,oc,sex\n4,0,male\n3,0,female\n2,0,male\n1,0,female\n")
            env = dict(os.environ, REPRO_LEGACY_DIR=str(root))
            result = subprocess.run([sys.executable, str(script), "4", "0", "--legacy"],
                                    env=env, capture_output=True, text=True)
            # Existing CLI has no explicit legacy flag; run that same impossible
            # instance without it to exercise the old successful exhaustion path.
            if "could not convert string to float: '--legacy'" in result.stderr:
                result = subprocess.run([sys.executable, str(script), "4", "0"],
                                        capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertNotIn("temps=", result.stdout, "reported a time for an invalid optimum")


if __name__ == "__main__":
    unittest.main()
