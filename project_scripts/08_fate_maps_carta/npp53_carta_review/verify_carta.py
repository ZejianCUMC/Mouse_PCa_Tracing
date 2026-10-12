import json, re, sys

PATH = '/mnt/7efffd2d-39a7-4d6a-87c7-773685593413/dropbox/Dropbox/share/Npp53_tracing/Npp53_tree_viewer_interactive.html'
html = open(PATH).read()

def extract_balanced(s, start_key):
    idx = s.find(start_key)
    if idx == -1:
        raise ValueError(f"key {start_key} not found")
    # find the opening brace/bracket right after the colon
    ci = s.find(':', idx) + 1
    while s[ci] in ' \t\n':
        ci += 1
    open_ch = s[ci]
    close_ch = {'{':'}', '[':']'}[open_ch]
    depth = 0
    in_str = False
    esc = False
    for j in range(ci, len(s)):
        c = s[j]
        if in_str:
            if esc:
                esc = False
            elif c == '\\':
                esc = True
            elif c == '"':
                in_str = False
            continue
        else:
            if c == '"':
                in_str = True
            elif c == open_ch:
                depth += 1
            elif c == close_ch:
                depth -= 1
                if depth == 0:
                    return s[ci:j+1]
    raise ValueError("unbalanced")

state_order_raw = extract_balanced(html, '"stateOrder"')
state_order = json.loads(state_order_raw)

palette_raw = extract_balanced(html, '"palettes"')
palettes = json.loads(palette_raw)
state_palette = palettes['state']

carta_raw = extract_balanced(html, '"carta"')
carta = json.loads(carta_raw)

print(f"Found {len(state_order)} states in stateOrder")
print(f"Found {len(state_palette)} states in palette")
print(f"Found {len(carta)} samples in carta")
print("Samples:", sorted(carta.keys()))

json.dump(carta, open('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review/carta_extracted.json','w'))
json.dump({'stateOrder': state_order, 'palette': state_palette}, open('/media/jrchristin/Backup10TB1/agent/1-project/12-Jingbo/2-Tracing/1-mouse_based/Npp53_tracing/analysis/evolution_features_v2/carta_fine12/biology_review/state_meta.json','w'))
