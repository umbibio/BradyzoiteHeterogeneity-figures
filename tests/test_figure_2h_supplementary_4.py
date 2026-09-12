"""Run: python -m unittest discover -s tests -v. No statistical reruns."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")  # Tests render without a desktop GUI (including on macOS).
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from bzfig import figure_2h_supplementary_4 as heatmaps
from bzfig.data import load_dataset


class RegistryTests(unittest.TestCase):
    def test_list_without_data(self):
        # --list is a registry question, so it must answer without a dataset.
        with tempfile.TemporaryDirectory() as empty:
            command = [sys.executable, str(ROOT / "scripts/make_figures.py"), "--list", "--data", empty]
            run = subprocess.run(command, capture_output=True, text=True, check=True)
        listed = run.stdout.splitlines()
        self.assertEqual(len(listed), 45)
        self.assertEqual(len(set(listed)), 45)
        for name in ("Figure_2H_correlation_heatmap", "Supplementary_4_CCC", "Supplementary_4_MCC"):
            self.assertIn(name, listed)


class OriginalHeatmapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.adata = load_dataset(ROOT / "data")
        cls.reference = json.loads((ROOT / "tests/reference_heatmaps.json").read_text())
        cls.phase = heatmaps.phase_matrices(cls.adata)
        cls.correlation = heatmaps.correlation_matrix(cls.adata)

    def test_original_phase_values_and_orders(self):
        for cycle in ("CCC", "MCC"):
            with self.subTest(cycle=cycle):
                expected = self.reference["supp4"][cycle]
                mean = np.array(expected["mean"])
                z = (mean-mean.mean(1, keepdims=True))/mean.std(1, ddof=1, keepdims=True)
                self.assertEqual(list(self.phase[cycle]["mean"].index), expected["gene_ids"])
                np.testing.assert_allclose(self.phase[cycle]["mean"], mean, atol=1e-6, rtol=0)
                np.testing.assert_allclose(self.phase[cycle]["z"], z, atol=1e-6, rtol=0)

    def test_original_correlations(self):
        expected = np.array(self.reference["fig2H"])
        np.testing.assert_allclose(self.correlation, expected, atol=1e-6, rtol=0)
        np.testing.assert_array_equal(np.round(self.correlation, 2), np.round(expected, 2))

    def test_no_mutation_of_shared_data(self):
        before = self.adata.layers["logcounts"].copy()
        heatmaps.normalized_expression(self.adata)
        self.assertEqual((before != self.adata.layers["logcounts"]).nnz, 0)

    def test_refuse_wrong_gene_universe(self):
        with self.assertRaisesRegex(ValueError, "8,170"):
            heatmaps.normalized_expression(self.adata[:, :100].copy())

    def test_h5ad_and_plain_inputs_agree(self):
        alternate = load_dataset(ROOT / "data", source="h5ad")
        np.testing.assert_allclose(heatmaps.correlation_matrix(alternate), self.correlation, atol=1e-12, rtol=0)

    def test_tables_are_written_from_results(self):
        with tempfile.TemporaryDirectory() as out:
            heatmaps.write_tables(Path(out), correlation=self.correlation, phase=self.phase)
            self.assertEqual(len(list(Path(out).glob("*.csv"))), 7)

    def test_original_manual_highlight_annotations(self):
        # Independent transcription of Supplementary information full.pdf p4.
        # The source has 19 CCC strips and 20 MCC strips, not 20 in both.
        import matplotlib.pyplot as plt
        from matplotlib.colors import to_hex
        shared = {"260250", "313040", "223050", "254910", "220440", "229200",
                  "249880", "219100", "261410", "216220", "240460", "281450",
                  "262730", "203710", "270330", "318470", "251740", "237090", "207900"}
        for cycle in ("CCC", "MCC"):
            with self.subTest(cycle=cycle):
                ids = shared | ({"315760"} if cycle == "MCC" else set())
                expected = {"TGME49_" + gene for gene in ids}
                self.assertEqual(set(heatmaps.METADATA["supp4_highlights"][cycle]), expected)
                fig = heatmaps.supplementary_4(self.phase[cycle], cycle)
                try:
                    ax = fig.axes[0]
                    patches = {p.get_gid().removeprefix("manual-highlight-"): p
                               for p in ax.patches if (p.get_gid() or "").startswith("manual-highlight-")}
                    self.assertEqual(set(patches), expected)
                    for gene, patch in patches.items():
                        self.assertEqual(to_hex(patch.get_facecolor()), "#ffcde6")
                        row = list(self.phase[cycle]["z"].index).index(gene)
                        self.assertEqual(patch.get_y(), row - 0.5)
                    self.assertTrue(all(t.get_color() == "black" for t in ax.texts))
                    np.testing.assert_array_equal(ax.images[0].get_array(), self.phase[cycle]["z"])
                finally:
                    plt.close(fig)


if __name__ == "__main__":
    unittest.main()
