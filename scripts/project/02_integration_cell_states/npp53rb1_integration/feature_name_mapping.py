"""Lossless stable-ID mapping to validated Seurat feature names."""

def seurat_feature_names(symbols):
    """R make.unique(..., sep='.') followed by Seurat underscore/pipe repair.

    Reserve original names before allocating suffixes, including genuine names
    ending in .1. Never merge distinct stable gene IDs with the same symbol.
    """
    reserved = set(symbols)
    seen = set()
    next_suffix = {}
    result = []
    for symbol in symbols:
        name = symbol
        if name in seen:
            suffix = next_suffix.get(symbol, 1)
            while f'{symbol}.{suffix}' in reserved:
                suffix += 1
            name = f'{symbol}.{suffix}'
            next_suffix[symbol] = suffix + 1
            reserved.add(name)
        seen.add(name)
        result.append(name.replace('_', '-').replace('|', '-'))
    if len(set(result)) != len(result):
        raise ValueError('Seurat character repair creates collisions; explicit mapping required')
    return result


def map_corrected_features(raw_ids, raw_symbols, exported_names,
                           corrected_ids, corrected_symbols):
    """Return existing unique names in corrected-matrix row order, via IDs."""
    if not (len(raw_ids) == len(raw_symbols) == len(exported_names)):
        raise ValueError('Raw reference and exported Seurat feature counts differ')
    if len(corrected_ids) != len(corrected_symbols):
        raise ValueError('Corrected feature IDs and symbols have different lengths')
    if len(set(raw_ids)) != len(raw_ids) or len(set(corrected_ids)) != len(corrected_ids):
        raise ValueError('Stable gene IDs must be unique')
    if seurat_feature_names(raw_symbols) != list(exported_names):
        raise ValueError('Raw reference does not reproduce exported Seurat feature names/order')
    lookup = dict(zip(raw_ids, zip(raw_symbols, exported_names)))
    result = []
    for gene_id, symbol in zip(corrected_ids, corrected_symbols):
        if gene_id not in lookup:
            raise ValueError(f'Corrected gene ID absent from original reference: {gene_id}')
        old_symbol, unique_name = lookup[gene_id]
        if old_symbol != symbol:
            raise ValueError(f'Conflicting symbols for stable gene ID {gene_id}')
        result.append(unique_name)
    if len(set(result)) != len(result):
        raise ValueError('Mapped features are not unique')
    return result
