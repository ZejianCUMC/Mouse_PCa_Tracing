"""Apply the user's cluster-based tumour definition; retain reporter calls separately."""
import csv
import json
from collections import Counter

HOST_CLUSTERS = {'4', '5', '6', '7'}
TREATMENTS = {'Intact': 0, 'Castration': 1, 'Castration+ENZA': 2}

def build_cluster_tumor_tables(review):
    with (review / 'marker_first_all_cell_annotations.csv').open(newline='') as f:
        source = list(csv.DictReader(f))
    cell_key = next(iter(source[0]))
    assert len(source) == 56686 and len({r[cell_key] for r in source}) == 56686
    counts = Counter()
    total_counts = Counter()
    records = []
    for row in source:
        cluster = row['primary_cluster']
        assert cluster in {str(i) for i in range(13)}
        treatment = row['treatment'].replace('Castration+Enzalutamide', 'Castration+ENZA')
        assert treatment in TREATMENTS
        key = (row['sample_id'], row['library'], treatment)
        is_tumor = cluster not in HOST_CLUSTERS
        counts[key] += int(is_tumor)  # Keep a zero row if a sample has no tumour cells.
        total_counts[key] += 1
        records.append(dict(cell_id=row[cell_key], sample_id=key[0], run_id=key[1],
                            treatment=treatment, harmony_cluster=cluster,
                            tumor_by_cluster=is_tumor,
                            reporter_positive=float(row['raw_GFP']) > 0 or float(row['raw_mCherry']) > 0))
    assert len(counts) == 22 and len({k[0] for k in counts}) == 22
    table = [dict(sample_id=k[0], run_id=k[1], treatment=k[2],
                  number_of_total_cells=total_counts[k], number_of_tumor_cells=n)
             for k, n in counts.items()]
    table.sort(key=lambda r: (TREATMENTS[r['treatment']], r['number_of_tumor_cells'], r['sample_id']))
    for name, rows in [('tumor_cells_by_sample.csv', table), ('cluster_tumor_cell_classification.csv', records)]:
        with (review / name).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    result = dict(definition='User-defined tumour classification: Harmony theta1 / Leiden 0.3 seed 19 cluster NOT IN 4,5,6,7; independent of reporter status',
                  excluded_clusters=sorted(HOST_CLUSTERS), tumor_clusters=[0,1,2,3,8,9,10,11,12],
                  tumor_cells=sum(counts.values()), non_tumor_cells=len(source)-sum(counts.values()),
                  tumor_reporter_positive=sum(r['tumor_by_cluster'] and r['reporter_positive'] for r in records),
                  tumor_reporter_negative=sum(r['tumor_by_cluster'] and not r['reporter_positive'] for r in records),
                  non_tumor_reporter_positive=sum(not r['tumor_by_cluster'] and r['reporter_positive'] for r in records),
                  samples=22, caveat='Operational cluster-based classification, not independent confirmation of malignancy. Cluster 9 remains included with its QC/stress flag.')
    assert result['tumor_cells'] == 43057 and result['non_tumor_cells'] == 13629
    assert result['tumor_reporter_positive'] == 39683 and result['tumor_reporter_negative'] == 3374
    assert result['non_tumor_reporter_positive'] == 129
    (review / 'cluster_tumor_summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
