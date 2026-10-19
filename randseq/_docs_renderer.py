"""How `show_doc` lays out an API entry in these docs: the summary, the call signature, then a
table of the parameters with their types, defaults and docments, and what the function returns.
Used by the docs build only (`renderer` in `pyproject.toml`); the package never imports it."""
import inspect, re, textwrap
from fastcore.docments import docments
from nbdev.doclinks import NbdevLookup
from nbdev.showdoc import BasicMarkdownRenderer

_empty = inspect.Parameter.empty

def _type(anno):
    "A type annotation as a reader expects it: `DataFrame`, `list | None`."
    if anno is _empty: return ""
    s = anno.__name__ if isinstance(anno, type) else str(anno)
    return re.sub(r"\b(?:\w+\.)+(\w+)", r"\1", s).replace("NoneType", "None")   # pandas.DataFrame -> DataFrame

def _cell(s, code=False):
    """Text for a pipe-table cell, on one line: plain text with pipes and angle brackets as HTML
    entities outside code spans, so they neither split the cell nor read as tags; or code (no
    default value here contains a pipe)."""
    s = str(s).replace("\n", " ")
    if code: return f"`{s}`" if s else s
    esc = lambda t: t.replace("|", "&#124;").replace("<", "&lt;").replace(">", "&gt;")
    # outside code spans only: inside one, an entity would print as written
    return "`".join(p if i % 2 else esc(p) for i, p in enumerate(s.split("`")))

def _signature(name, sym):
    "`name(a, b=1, ...)`, wrapped to the page."
    try: params = inspect.signature(sym).parameters.values()
    except (TypeError, ValueError): return f"{name}(...)"
    parts = []
    for p in params:
        if p.name == "self": continue
        star = {p.VAR_POSITIONAL: "*", p.VAR_KEYWORD: "**"}.get(p.kind, "")
        parts.append(star + p.name + ("" if p.default is _empty else f"={p.default!r}"))
    return textwrap.fill(f"{name}({', '.join(parts)})", width=80, subsequent_indent=" " * (len(name) + 1),
                         break_long_words=False, break_on_hyphens=False)

class DocsRenderer(BasicMarkdownRenderer):
    "An nbdev `show_doc` renderer: summary, signature, parameter table."
    def _repr_markdown_(self):
        md = "---\n\n"
        src = NbdevLookup().code(self.fn)
        if src: md += f'[source]({src}){{target="_blank" style="float:right; font-size:smaller"}}\n\n'
        md += f"{'#' * self.title_level} {self.nm}\n\n"
        if self.docs: md += f"{self.docs}\n\n"
        if not callable(self.sym): return md
        md += f"```python\n{_signature(self.nm.split('.')[-1], self.sym)}\n```\n\n"
        try: dm = docments(self.sym, full=True)
        except Exception: return md
        rows = [(n, d) for n, d in dm.items() if n not in ("self", "return")]
        ret = dm.get("return", {})
        if rows or ret.get("docment"):
            md += "| Parameter | Type | Default | Description |\n|---|---|---|---|\n"
            for n, d in rows:
                default = "" if d["default"] is _empty else repr(d["default"])
                md += (f"| `{n}` | {_cell(_type(d['anno']))} | {_cell(default, code=True)} "
                       f"| {_cell(d['docment'] or '')} |\n")
            if ret.get("docment") or ret.get("anno", _empty) is not _empty:
                md += f"| *returns* | {_cell(_type(ret.get('anno', _empty)))} | | {_cell(ret.get('docment') or '')} |\n"
        return md
