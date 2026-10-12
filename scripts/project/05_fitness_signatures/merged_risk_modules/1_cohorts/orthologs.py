"""opt5 shared module — mouse -> human ortholog mapping (audit fixes 2026-10-09).
 1. MGI HOM_MouseHumanSequence homology classes: a mouse symbol maps to ALL human symbols of its class. opt3/opt4 kept only the first
    one (drop_duplicates), which was arbitrary for 434 tested one-to-many genes.
 2. Rescue of outdated mouse symbols. The 10x reference behind shared_counts.h5ad uses older MGI symbols (Atp5a1, Atp5b, Aars, Acpp,
    Arntl, ...), which are absent from the current HOM file, so about 2.5% of tested genes (including the ATP-synthase subunits) had
    no human ortholog. The uppercase symbol is looked up in HGNC (a current symbol, else a unique previous/alias symbol) and accepted
    only if that human gene belongs to a mouse-human HOM class (e.g. Atp5a1 -> ATP5A1 -> ATP5F1A). Rescues are reported by rescues().
Tracing transgenes (mCherry, Cas9, GFP, Cre) have no ortholog and are excluded upstream (02 / 03)."""
import pandas as pd
HOM = "/media/jrchristin/Backup10TB1/agent/1-project/11-prostate_atlas/22Science_Tang_CRPC_ATAC/resources/orthology/HOM_MouseHumanSequence.rpt"
HGNC = "/media/jrchristin/Backup10TB1/agent/1-project/25-CGC/data/references/hgnc/hgnc_complete_set_20260930.txt"
TRANSGENES = {"mCherry", "Cas9", "GFP", "Cre"}
_C = {}
def _load():
    if _C: return _C
    hom = pd.read_csv(HOM, sep="\t")
    mo = hom[hom["Common Organism Name"] == "mouse, laboratory"][["DB Class Key", "Symbol"]]
    hu = hom[hom["Common Organism Name"] == "human"][["DB Class Key", "Symbol"]]
    mh = mo.merge(hu, on="DB Class Key", suffixes=("_m", "_h"))
    _C["m2h"] = mh.groupby("Symbol_m").Symbol_h.apply(lambda s: sorted(set(s))).to_dict(); _C["hom_h"] = set(hu.Symbol)
    hg = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False, usecols=["symbol", "prev_symbol", "alias_symbol"])
    cur = set(hg.symbol); alias = {}
    for s, p, a in zip(hg.symbol, hg.prev_symbol.fillna(""), hg.alias_symbol.fillna("")):
        for x in set(p.replace('"', "").split("|") + a.replace('"', "").split("|")):
            if x and x not in cur: alias.setdefault(x, set()).add(s)
    _C["cur"] = cur; _C["a2c"] = {k: next(iter(v)) for k, v in alias.items() if len(v) == 1}
    return _C
def mouse_to_human(g):
    """(list of human symbols, how) for one mouse symbol; how = HOM | HGNC-rescue | none."""
    c = _load()
    if g in c["m2h"]: return c["m2h"][g], "HOM"
    if g in TRANSGENES: return [], "none"
    u = g.upper(); s = u if u in c["cur"] else c["a2c"].get(u)
    if s and s in c["hom_h"]: return [s], "HGNC-rescue"
    return [], "none"
def to_human(genes):
    return sorted({h for g in genes for h in mouse_to_human(g)[0]})
def rescues(genes):
    rows = [(g, *mouse_to_human(g)) for g in sorted(set(genes))]
    return pd.DataFrame([(g, ";".join(h), how) for g, h, how in rows], columns=["mouse", "human", "how"])
