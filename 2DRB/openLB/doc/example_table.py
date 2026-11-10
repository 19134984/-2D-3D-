#!/usr/bin/env python3
#
# Pre-render helper for the OpenLB user guide.
#
# Builds the tabular overview of the examples -- one row per example, one
# column per feature -- in two forms: as HTML for the Examples page of the
# website (examples/index_web.qmd), and as LaTeX for the A3 pages that open
# the examples chapter of the PDF (examples/index_pdf.qmd).
#
# The table is assembled from two sources:
#
#   * examples/_includes.qmd, written by link_examples.sh, which must therefore
#     run first. It lists every documented example in the order of the user
#     guide, grouped by category, so rows, categories and the links to the
#     example pages all follow the source tree without further maintenance.
#   * the front matter of each example's own .qmd file, which says which
#     features the example covers, using the keys of FEATURES below:
#
#         ---
#         features: [transient, benchmark, primitives]
#         recommended-for-beginners: [primitives]
#         ---
#
#     Features listed under recommended-for-beginners are shown in green and
#     must also appear under features. An example without a features key is
#     shown as not classified yet ("?"); one that covers none of them says
#     features: []. Quarto drops the front matter of included files, so it
#     does not show up in the PDF, and on the website it is plain metadata.
#
# In the PDF, the folder column refers to the section of each category and the
# example names link to the section of each example. The targets are the
# labels excat:<category> and ex:<path> that link_examples.sh sets up, the
# latter with the help of include/example_labels.lua.
#
# Both outputs, examples/_table.qmd and examples/_table_pdf.qmd, are generated
# and listed in .gitignore. The script is registered as a pre-render script in
# _quarto.yml and _quarto-book.yml after link_examples.sh, so a plain
# `quarto render` regenerates them automatically. It is idempotent and can
# also be run by hand.

import html
import re
import sys
from pathlib import Path

doc_dir = Path(__file__).resolve().parent
examples_dir = doc_dir / "examples"
includes = examples_dir / "_includes.qmd"
web_table = examples_dir / "_table.qmd"
pdf_table = examples_dir / "_table_pdf.qmd"

# The PDF table is split over as many A3 pages as it takes, at most this many
# rows each, and never within a category. A3 landscape fits a few more at the
# font size used below; the margin keeps the legend on the page.
PDF_ROWS_PER_PAGE = 50

# The columns of the table, in order: the key used in the front matter of the
# examples and the column heading. Add an entry here to add a column.
FEATURES = [
    ("turbulent", "turbulent"),
    ("thermal", "thermal"),
    ("multiComponent", "multi Component"),
    ("multiPhase", "multi Phase"),
    ("particles", "particles"),
    ("porousMedia", "porous Media"),
    ("transient", "transient flow"),
    ("benchmark", "benchmark"),
    ("showcase", "showcase"),
    ("stl", "STL geometry"),
    ("primitives", "geometry primitives"),
    ("checkpointing", "check-pointing"),
]


def log(msg):
    print(f"example_table: {msg}", file=sys.stderr)


def read_examples():
    """Return [(category, [example path, ...]), ...] in user guide order."""
    if not includes.is_file():
        log(f"no {includes.name} -- run link_examples.sh first")
        sys.exit(1)
    categories = []
    for line in includes.read_text().splitlines():
        if line.startswith("## "):
            categories.append((line[3:].strip(), []))
            continue
        m = re.match(r"\{\{< include (\S+)\.qmd >\}\}", line)
        # Hand written category introductions (_intro.qmd) are not examples.
        if m and not Path(m.group(1)).name.startswith("_"):
            categories[-1][1].append(m.group(1))
    return [(c, exs) for c, exs in categories if exs]


def front_matter_lists(path):
    """Return the list valued keys of the YAML front matter of a .qmd file.

    Only the two forms used for the feature lists are understood -- a flow
    list, key: [a, b], and a block list of "- a" lines -- which keeps the
    script free of a YAML library. Other keys are skipped.
    """
    lines = path.read_text().splitlines()
    if not lines or lines[0].strip() != "---":
        return {}
    lists, key = {}, None
    for line in lines[1:]:
        if line.strip() in ("---", "..."):
            break
        if re.match(r"\s*(#.*)?$", line):  # blank or comment only
            continue
        item =re.match(r"\s+-\s*(.*?)\s*(#.*)?$", line)
        if item and key:
            lists[key].append(item.group(1).strip("'\""))
            continue
        m = re.match(r"([\w-]+):\s*(.*?)\s*(#.*)?$", line)
        key = None
        if not m:
            continue
        name, value = m.group(1), m.group(2)
        if value.startswith("[") and value.endswith("]"):
            lists[name] = [v.strip().strip("'\"")
                           for v in value[1:-1].split(",") if v.strip()]
        elif value == "":
            key = name
            lists[name] = []
    return lists


def read_features(examples):
    """Return {example: [mark, ...]} with marks "", "x", "b" or "?"."""
    keys = [k for k, _ in FEATURES]
    features = {}
    todo = []
    for _, exs in examples:
        for ex in exs:
            meta = front_matter_lists(examples_dir / f"{ex}.qmd")
            if "features" not in meta:
                todo.append(ex)
                features[ex] = ["?"] * len(keys)
                continue
            have = set(meta["features"])
            beginner = set(meta.get("recommended-for-beginners", []))
            for k in sorted((have | beginner) - set(keys)):
                log(f"{ex}: unknown feature '{k}' (known: {', '.join(keys)})")
            for k in sorted(beginner - have):
                log(f"{ex}: '{k}' is recommended for beginners but not "
                    "listed under features")
            features[ex] = ["b" if k in beginner else "x" if k in have else ""
                            for k in keys]
    if todo:
        log(f"{len(todo)} example(s) without a features list in their front "
            "matter: " + ", ".join(todo))
    return features


def split_path(example):
    """laminar/cavity2d/cavity2d -> ("", "cavity2d"); keeps intermediate
    directories. Every page sits in a directory of its own, named after the
    example, which is left out."""
    parts = example.split("/")
    prefix = "/".join(parts[1:-2])
    return (prefix + "/" if prefix else ""), parts[-1]


# ---------------------------------------------------------------- website --

def web_cell(mark, dark):
    shade = "dark" if dark else ""
    cls = {
        "": "bg-darkgray" if dark else "bg-lightgray",
        "x": f"bg-royalblue{shade}",
        "b": "bg-green",
        "?": "bg-unknown",
    }[mark]
    title = {"b": ' title="recommended for beginners"',
             "?": ' title="not classified yet"'}.get(mark, "")
    text = "?" if mark == "?" else "&nbsp;"
    return f'<td class="{cls}"{title}>{text}</td>'


def write_web(examples, features):
    e = html.escape
    out = []
    out.append("<!-- Generated by example_table.py -- do not edit. -->")
    out.append("```{=html}")
    # The legend sits above the table, outside the horizontally scrolling
    # container, so that it stays in view.
    out.append('<table class="legend-table footnotesize">')
    out.append("<tr>")
    out.append('<td class="bg-royalblue">&nbsp;</td>'
               '<td class="bg-royalbluedark">&nbsp;</td>'
               "<td>example includes relevant subject</td>")
    out.append('<td class="bg-green">&nbsp;</td>'
               "<td>example includes relevant subject and is recommended "
               "for beginning</td>")
    out.append('<td class="bg-unknown">?</td>'
               "<td>not classified yet</td>")
    out.append("</tr>")
    out.append("</table>")
    out.append('<div class="table-wrap footnotesize">')
    out.append('<table class="matrix">')
    out.append("<thead>")
    out.append('<tr class="thick-bottom">')
    out.append("<th>folder</th>")
    out.append('<th class="left" style="min-width: 240px;">example</th>')
    for _, label in FEATURES:
        out.append(f"<th>{e(label)}</th>")
    out.append("</tr>")
    out.append("</thead>")
    out.append("<tbody>")
    for category, exs in examples:
        out.append(f"<!-- {e(category)} ({len(exs)}) -->")
        for i, ex in enumerate(exs):
            dark = i % 2 == 0
            bg = "bg-darkgray" if dark else "bg-lightgray"
            out.append('<tr class="thick-top">' if i == 0 else "<tr>")
            if i == 0:
                out.append(f'<td class="folder" rowspan="{len(exs)}">'
                           f"{e(category)}</td>")
            prefix, name = split_path(ex)
            prefix = f'<span class="ex-path">{e(prefix)}</span>' if prefix else ""
            # The table sits on examples/index_web.qmd, next to the mirrored
            # example pages, so the link is relative to examples/.
            out.append(f'<td class="left {bg}"><a href="{e(ex)}.html">'
                       f"{prefix}{e(name)}</a></td>")
            for mark in features[ex]:
                out.append(web_cell(mark, dark))
            out.append("</tr>")
    out.append("</tbody>")
    out.append("</table>")
    out.append("</div>")
    out.append("```")
    web_table.write_text("\n".join(out) + "\n")

# -------------------------------------------------------------------- PDF --

def tex(text):
    """Escapes the characters LaTeX treats specially."""
    special = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
               "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
               "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(special.get(c, c) for c in text)


def tex_head(label):
    """A column heading over two lines, like "multi / Component"."""
    if " " in label:
        first, rest = label.split(" ", 1)
    elif "-" in label:
        first, rest = label.split("-", 1)
        first += "-"
    else:
        # One line, vertically centred against the two line headings.
        return r"\begin{tabular}[c]{@{}c@{}}" + tex(label) + r"\end{tabular}"
    return (r"\begin{tabular}[c]{@{}c@{}}" + tex(first) + r"\\" + tex(rest)
            + r"\end{tabular}")


def paginate(examples):
    """Splits the categories over as few pages as possible, evenly.

    The categories stay in order and whole. Among the splits into the fewest
    pages of at most PDF_ROWS_PER_PAGE rows, the one with the shortest
    longest page is taken, so that the pages come out about equally full.
    """
    sizes = [len(exs) for _, exs in examples]
    for c, n in zip(examples, sizes):
        if n > PDF_ROWS_PER_PAGE:
            log(f"category {c[0]} has {n} examples, more than fit on a page "
                f"({PDF_ROWS_PER_PAGE}); its page will overflow")
    limit = max([PDF_ROWS_PER_PAGE] + sizes)
    total = sum(sizes)
    pages = -(-total // limit)
    while True:
        # best[k][i]: smallest possible longest page when the first i
        # categories go onto k pages, and where the last of them starts.
        inf = float("inf")
        best = [[(inf, 0)] * (len(sizes) + 1) for _ in range(pages + 1)]
        best[0][0] = (0, 0)
        for k in range(1, pages + 1):
            for i in range(1, len(sizes) + 1):
                rows = 0
                for j in range(i, 0, -1):
                    rows += sizes[j - 1]
                    if rows > limit:
                        break
                    cand = max(best[k - 1][j - 1][0], rows)
                    if cand < best[k][i][0]:
                        best[k][i] = (cand, j - 1)
        if best[pages][len(sizes)][0] < inf:
            break
        pages += 1
    cuts, i = [], len(sizes)
    for k in range(pages, 0, -1):
        j = best[k][i][1]
        cuts.append((j, i))
        i = j
    return [examples[a:b] for a, b in reversed(cuts)]


def pdf_cell(mark, dark):
    color = {
        "": "darkgray" if dark else "lightgray",
        "x": "RoyalBluedark" if dark else "RoyalBlue",
        "b": "olbBeginner",
        "?": "olbUnknown",
    }[mark]
    return rf"\cellcolor{{{color}}}" + ("?" if mark == "?" else "")


def pdf_page(examples, features):
    cols = len(FEATURES)
    out = []
    out.append(r"\begin{tabularx}{\textwidth}{|l|l|*{%d}{>{\centering\arraybackslash}X|}}" % cols)
    out.append(r"\noalign{\hrule height 1pt}")
    heads = [tex_head("folder (section)"), tex_head("example")]
    heads += [tex_head(label) for _, label in FEATURES]
    out.append(" & ".join(heads) + r" \\")
    out.append(r"\noalign{\hrule height 1pt}")
    for category, exs in examples:
        for i, ex in enumerate(exs):
            dark = i % 2 == 0
            folder = ""
            if i == 0:
                # Centred over the rows of the category, as \multirow would.
                folder = (rf"\olbFolder{{{len(exs)}}}{{{tex(category)} "
                          rf"(\ref{{excat:{category}}})}}")
            prefix, name = split_path(ex)
            if prefix:
                # Nested examples show the directories in between, which may
                # break the line after a slash.
                prefix = r"/\allowbreak ".join(tex(p) for p in prefix.split("/"))
                prefix = rf"{{\color{{gray}}{prefix}}}\allowbreak "
            example = (rf"\cellcolor{{{'darkgray' if dark else 'lightgray'}}}"
                       rf"\hyperref[ex:{ex}]{{{prefix}{tex(name)}}}")
            cells = [folder, example] + [pdf_cell(m, dark) for m in features[ex]]
            out.append(" & ".join(cells) + r" \\")
        out.append(r"\noalign{\hrule height 1pt}")
    out.append(r"\end{tabularx}")
    out.append(r"\par\vspace{2mm}\noindent")
    out.append(r"\olbKey{RoyalBlue}\,\olbKey{RoyalBluedark}~example includes relevant subject\qquad")
    out.append(r"\olbKey{olbBeginner}~example includes relevant subject and is recommended for beginning\qquad")
    out.append(r"\olbKey{olbUnknown}~not classified yet")
    return out


def write_pdf(examples, features):
    pages = paginate(examples)
    out = []
    out.append("<!-- Generated by example_table.py -- do not edit. -->")
    out.append("```{=latex}")
    out.append(r"% The A3 page layout around this table is set up by index_pdf.qmd.")
    out.append(r"\begingroup")
    out.append(r"\definecolor{olbBeginner}{RGB}{0,160,0}")
    out.append(r"\definecolor{olbUnknown}{RGB}{255,233,168}")
    out.append(r"\hypersetup{linkcolor=black}")
    out.append(r"\singlespacing\footnotesize")
    out.append(r"\renewcommand{\arraystretch}{1.25}")
    out.append(r"\newlength{\olbRow}")
    out.append(r"\setlength{\olbRow}{\dimexpr\ht\strutbox+\dp\strutbox\relax}")
    out.append(r"\setlength{\olbRow}{1.25\olbRow}")
    # #1 rows high: lowered from the first row by half the rows below it.
    out.append(r"\newcommand{\olbFolder}[2]{\smash{\raisebox{-0.5\dimexpr\olbRow*(#1-1)\relax}{#2}}}")
    out.append(r"\newcommand{\olbKey}[1]{\textcolor{#1}{\rule[-0.2em]{1.2em}{1em}}}")
    for n, page in enumerate(pages):
        if n:
            out.append(r"\clearpage")
        out.append(r"\thispagestyle{empty}")
        out.append(r"\noindent")
        out.extend(pdf_page(page, features))
    out.append(r"\clearpage")
    out.append(r"\endgroup")
    out.append("```")
    pdf_table.write_text("\n".join(out) + "\n")
    return pages


def main():
    examples = read_examples()
    features = read_features(examples)
    write_web(examples, features)
    pages = write_pdf(examples, features)
    rows = sum(len(exs) for _, exs in examples)
    print(f"example_table: wrote {web_table.name} and {pdf_table.name} "
          f"({rows} examples, {len(pages)} A3 pages: "
          + ", ".join(str(sum(len(e) for _, e in p)) for p in pages)
          + " rows)")


if __name__ == "__main__":
    main()
