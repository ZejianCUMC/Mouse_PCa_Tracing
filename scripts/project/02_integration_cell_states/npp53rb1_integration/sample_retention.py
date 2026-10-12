"""User-specified per-mouse cutoff, applied after RNA QC and before integration."""
MIN_VALID_CELLS = 400
TREATMENT_ORDER = {'Intact': 0, 'Castration': 1, 'Castration+ENZA': 2,
                   'Castration+Enzalutamide': 2}

def retain_sample(count, minimum=MIN_VALID_CELLS):
    return int(count) >= minimum

def sample_sort_key(row):
    return (TREATMENT_ORDER[row['treatment']], int(row['cells']), row['sample_id'])
