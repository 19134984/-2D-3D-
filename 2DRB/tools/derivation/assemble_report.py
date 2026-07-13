"""Assemble the complete Markdown report from reviewed source chapters."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "docs/lbm-cde-trt-derivation.md"
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


def assemble() -> str:
    """Return the normalized, deterministic complete report text."""
    return "\n\n".join(
        part.read_text(encoding="utf-8").rstrip() for part in PARTS
    ) + "\n"


def main() -> None:
    OUTPUT.write_text(assemble(), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
