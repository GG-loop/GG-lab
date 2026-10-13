#!/usr/bin/env python3
"""Copy export_snippet.js into every viewer page source between its markers."""
from pathlib import Path

HERE = Path(__file__).parent
START, END = "// >>> export_snippet.js", "// <<< export_snippet.js"
snippet = (HERE / "export_snippet.js").read_text().rstrip("\n")
for name in ["organic.html", "template.html", "feasibility_template.html", "sculpture_template.html"]:
    p = HERE / name
    s = p.read_text()
    a, b = s.index(START), s.index(END)
    s = s[:a] + START + "\n" + snippet + "\n" + s[b:]
    p.write_text(s)
    print("updated", name)
