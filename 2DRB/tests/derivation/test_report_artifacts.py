from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "docs/lbm-cde-trt-derivation.md"
PARTS = (
    ROOT / "docs/derivation/chapters/01-overview.md",
    ROOT / "docs/derivation/chapters/02-trt-source-design.md",
    ROOT / "docs/derivation/chapters/03-effective-rates.md",
    ROOT / "docs/derivation/chapters/04-second-order-recovery.md",
    ROOT / "docs/derivation/chapters/05-boundary-magic.md",
    ROOT / "docs/derivation/chapters/06-d2q9-fourth-order.md",
    ROOT / "docs/derivation/chapters/07-parameter-compatibility.md",
    ROOT / "docs/derivation/chapters/08-algorithm.md",
    ROOT / "docs/derivation/chapters/A-d2q5-verifier.md",
    ROOT / "docs/derivation/review-matrix.md",
    ROOT / "docs/derivation/references.md",
)


class CompleteMarkdownReportTests(unittest.TestCase):
    def test_primary_markdown_is_the_exact_full_report_assembly(self):
        expected = "\n\n".join(
            part.read_text(encoding="utf-8").rstrip() for part in PARTS
        ) + "\n"
        self.assertEqual(REPORT.read_text(encoding="utf-8"), expected)

    def test_primary_markdown_contains_detailed_derivation_landmarks(self):
        report = REPORT.read_text(encoding="utf-8")
        for landmark in (
            "逐分量 BGK 极限",
            "一般 ABB 独立 jets",
            "D2Q9 LBM-CDE-TRT 冻结系数四阶等效方程",
            "参数化 D2Q5 四阶参考验证器",
            "需求—证据—验证矩阵",
        ):
            with self.subTest(landmark=landmark):
                self.assertIn(landmark, report)

    def test_build_entrypoint_enforces_repeatable_metadata_and_assembly(self):
        build = (ROOT / "tools/derivation/build_report.ps1").read_text(
            encoding="utf-8"
        )
        tex = (ROOT / "docs/lbm-cde-trt-derivation.tex").read_text(
            encoding="utf-8"
        )
        self.assertIn("assemble_report.py", build)
        self.assertIn("SOURCE_DATE_EPOCH", build)
        self.assertIn("FORCE_SOURCE_DATE", build)
        self.assertNotIn(r"\today", tex)


if __name__ == "__main__":
    unittest.main()
