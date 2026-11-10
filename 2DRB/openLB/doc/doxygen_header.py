#!/usr/bin/env python3
"""Build a doxygen HTML header carrying the rendered Quarto navbar.

The navbar markup and its links are taken from the rendered user guide rather
than written out again here, so the API reference keeps the same top bar as the
guide it is published next to. Paths are rewritten one level down, because the
API reference is served from <site>/api/ while the guide sits at <site>/.
"""
import re
import subprocess
import tempfile
import sys
from pathlib import Path

doc = Path(__file__).resolve().parent
site_index = doc / "ug" / "index.html"
out = doc / "_doxygen_header.html"

if not site_index.exists():
    sys.exit(f"doxygen_header: no rendered guide at {site_index}")

html = site_index.read_text(encoding="utf8", errors="replace")

# The whole header element, not just the <nav> inside it: Quarto's stylesheet
# keys the fixed positioning and spacing off this wrapper and off the body
# classes copied below, so lifting the <nav> alone renders a different bar.
start = html.find('<header id="quarto-header"')
end = html.find("</header>", start)
if start < 0 or end < 0:
    sys.exit("doxygen_header: no quarto header found in the rendered guide")
nav = html[start:end + len("</header>")]

body_class = re.search(r'<body([^>]*)class="([^"]*)"', html)
body_class = body_class.group(2) if body_class else "nav-fixed quarto-light"

# The search box indexes the guide, not the API reference, so it would lie here.
nav = re.sub(r'<div id="quarto-search".*?</div>\s*', "", nav, flags=re.S)

# The guide's own stylesheets, minus the ones that only matter for its content.
css = [c for c in re.findall(r'<link[^>]+href="([^"]+\.css)"', html)
       if "pseudocode" not in c and "glightbox" not in c and "tippy" not in c]

def down(path: str) -> str:
    """Rewrite a link from the guide root to one level below it."""
    return "../" + path[2:] if path.startswith("./") else (
        "../" + path if not path.startswith(("http", "#", "/")) else path)

nav = re.sub(r'(href|src)="(\./[^"]*)"', lambda m: f'{m.group(1)}="{down(m.group(2))}"', nav)

# Start from doxygen's own header so its scripts and placeholders stay intact.
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run(["doxygen", "-w", "html", f"{tmp}/h.html", f"{tmp}/f.html",
                    f"{tmp}/s.css"], check=True, capture_output=True)
    tpl = Path(tmp, "h.html").read_text(encoding="utf8")

# doxygen-awesome switches to dark on prefers-color-scheme unless the document
# opts out with .light-mode. The guide has no dark variant (Quarto renders it
# quarto-light), so a dark-preference visitor would otherwise get dark chrome
# around a white content pane. Pin the reference to light as well.
tpl = re.sub(r"<html\b", '<html class="light-mode"', tpl, count=1)

links = "\n".join(f'<link href="{down(c)}" rel="stylesheet">' for c in css)

# Without these the bar is inert markup: the hamburger and the Documentation
# dropdown are Bootstrap behaviours, and the pinning is Quarto's. Search is
# deliberately left out along with its box.
wanted = ("bootstrap/bootstrap.min.js", "quarto-nav/headroom.min.js",
          "quarto-nav/quarto-nav.js")
scripts = [j for j in re.findall(r'<script[^>]+src="([^"]+\.js)"', html)
           if any(j.endswith(w) for w in wanted)]
links += "\n" + "\n".join(
    f'<script src="{down(j)}" defer></script>' for j in scripts)
tpl = tpl.replace("</head>", f"{links}\n</head>")
tpl = re.sub(r"<body[^>]*>", f'<body class="{body_class}">\\n' + nav.replace("\\\\", "\\\\\\\\"),
             tpl, count=1)

out.write_text(tpl, encoding="utf8")
print(f"doxygen_header: wrote {out.relative_to(doc.parent)} ({len(tpl)} bytes)")
