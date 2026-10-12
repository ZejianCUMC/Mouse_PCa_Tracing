"""Readable MSigDB Hallmark names: HALLMARK_KRAS_SIGNALING_DN -> 'Hallmark KRAS signaling (down)', keeping gene / pathway acronyms."""
KEEP = {"KRAS", "MYC", "E2F", "G2M", "MTORC1", "TNFA", "NFKB", "IL2", "STAT5", "IL6", "JAK", "STAT3", "TGF", "UV", "P53", "PI3K", "AKT", "MTOR", "WNT", "DNA", "UPR", "ROS"}
PRETTY = {"MTORC1": "mTORC1", "MTOR": "mTOR", "TNFA": "TNFα", "NFKB": "NF-κB", "V1": "v1", "V2": "v2", "DN": "(down)", "UP": "(up)"}
def hallmark(x):
    w = x.replace("HALLMARK_", "").split("_")
    return "Hallmark " + " ".join(PRETTY.get(t, t if t in KEEP else t.lower()) for t in w)
