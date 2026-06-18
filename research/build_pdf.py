"""Build a submission PDF from MANUSCRIPT.md: markdown -> styled HTML.

Writes research/MANUSCRIPT.html; the PDF is then produced by the weasyprint CLI
(its Python import is broken in this environment, the binary works):

    python3 research/build_pdf.py
    weasyprint research/MANUSCRIPT.html research/MANUSCRIPT.pdf

A serif academic layout with a wide-coverage font fallback so the unicode maths
(ᵀ, ℝ, μ, λ, ≤, ∈, Σ, ², 𝟙) renders. Figure 1 (fig_scaling.png) is embedded.
"""
import markdown
from pathlib import Path

R = Path("research")
md = (R / "MANUSCRIPT.md").read_text()

# Embed Figure 1 above its caption (the .md only references the PDF path).
md = md.replace(
    "## Figure\n\n**Figure 1.**",
    "## Figure\n\n![Figure 1](fig_scaling.png)\n\n**Figure 1.**",
)

body = markdown.markdown(md, extensions=["tables", "sane_lists", "attr_list"])

CSS = """
@page { size: A4; margin: 2cm 2cm 2.2cm 2cm;
        @bottom-center { content: counter(page); font: 9pt Georgia, serif; color:#666; } }
body { font-family: 'Times New Roman', Georgia, 'DejaVu Serif', 'DejaVu Sans', serif;
       font-size: 10.5pt; line-height: 1.42; color:#111; text-align: justify; }
h1 { font-size: 16pt; text-align:center; line-height:1.25; margin:0 0 0.15em; }
h1 + p { text-align:center; margin:0.1em 0; }           /* author / affiliation */
h2 { font-size: 12.5pt; border-bottom:1px solid #ccc; padding-bottom:2px;
     margin:1.25em 0 0.4em; page-break-after: avoid; }
h3 { font-size: 11pt; margin:1em 0 0.3em; page-break-after: avoid; }
p { margin:0.45em 0; }
table { border-collapse: collapse; width:100%; font-size:9.3pt; margin:0.7em 0;
        page-break-inside: avoid; }
th, td { border:1px solid #bbb; padding:3px 6px; text-align:left; }
th { background:#f0f0f0; }
code { font-family:'DejaVu Sans Mono', monospace; font-size:9pt; background:#f5f5f5;
       padding:0 2px; }
img { max-width:100%; display:block; margin:0.6em auto; }
a { color:#1a4e8a; text-decoration:none; }
blockquote { color:#444; border-left:3px solid #ddd; margin:0.6em 0; padding:0.1em 0.8em; }
"""

html = (
    "<!doctype html><html><head><meta charset='utf-8'>"
    f"<style>{CSS}</style></head><body>{body}</body></html>"
)
(R / "MANUSCRIPT.html").write_text(html)
print(f"wrote research/MANUSCRIPT.html ({len(body)} chars of body)")
