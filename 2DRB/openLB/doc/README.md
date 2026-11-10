# OpenLB documentation sources

The rendered documentation lives at <https://openlb.gitlab.io/olb/> — the user
guide at the root and the API reference under `/api`. It is republished by the
`pages` job in `.gitlab-ci.yml`, and access is restricted to project members.

Read it there — GitLab renders only a fixed set of markdown extensions, so the
`.qmd` files below are shown as source.

## Layout

| Path | Contents |
| --- | --- |
| `index.qmd`, `about.qmd`, `cite.qmd`, … | Landing and standalone pages |
| `content/` | The user guide chapters |
| `examples/` | One page per example, plus the generated index |
| `examples/example.qmd` | Template and guide for writing the page of a new example |
| `_quarto.yml` | Website profile (the default) |
| `_quarto-book.yml` | PDF book profile |
| `fig/`, `include/`, `styles.css` | Assets, header partials, theme overrides |
| `doxygen_header.py`, `doxygen-quarto.css` | Give the API reference this site's navbar |

## Building

Everything runs through the `docs` environment declared in the flake, which
pins Quarto, a TeX distribution, Doxygen and the helper tools:

    nix develop .#docs --command make userguide   # HTML site into doc/ug
    nix develop .#docs --command make docs        # the above plus the API reference

Both are run from the repository root. Inside `doc/` the individual targets are:

    make          # the PDF, then the website
    make pdf      # doc/ug_pdf/OpenLB-User-Guide.pdf
    make website  # doc/ug/index.html
    make clean    # removes both renders and everything they generate

Open `doc/ug/index.html` in a browser to read the result.

Quarto applies one profile per invocation, so the PDF and the website cannot be
rendered in a single call. The order matters: the website only shows the PDF
download button next to the GitLab icon in its navbar if the PDF is already on
disk, so render the book first. Plain `make` does this for you.

### Without Nix

Install Quarto (1.8.26 or newer), then `quarto install tinytex` for the PDF and
`rsvg-convert` for figure conversion (`apt install librsvg2-bin` on Debian and
Ubuntu). Run `quarto check` and confirm it picked up the TeX installation you
intend. The makefile targets above work the same way.

## Generated files

Several parts of the tree are produced at render time and must not be edited by
hand or committed:

- `link_examples.sh` mirrors each example's documentation into a directory of
  its own below `doc/examples/`, as symlinks to the page and, if the example has
  one, to its `doc/` directory of figures — the pages themselves live next to
  their code under `examples/<category>/` — and writes the
  `examples/_includes.qmd` list and the generated section of the sidebar in
  `_quarto.yml`.
- `example_table.py` builds the overview table of the examples, as HTML for
  the Examples page and as LaTeX for the A3 pages at the start of the examples
  chapter of the PDF, from the example list written by `link_examples.sh` and
  the front matter of each example's `.qmd`, which lists the features it
  covers:

      ---
      features: [transient, benchmark, primitives]
      recommended-for-beginners: [primitives]
      ---

  The feature keys and column headings are defined at the top of the script.
  An example without `features` is shown with `?` in every column.

  In the PDF, the table refers to the sections of the categories and links to
  those of the examples through the labels `excat:<category>` and
  `ex:<path>`. `link_examples.sh` writes them into the include list, and
  `include/example_labels.lua` places each example's label after its heading.
- `ug_pdf_download.sh` copies the finished PDF into the site and writes the
  navbar download button.

These scripts also run as Quarto pre-render steps, so a plain `quarto
render` regenerates them too.

The `_extensions/` directory vendors the `quarto-pseudocode` extension, which
renders the algorithm listings. It is checked in; run
`quarto add leovan/quarto-pseudocode` only to update it.
