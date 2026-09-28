#!/bin/bash
# =============================================================================
# hashtags.sh — writes the two small metadata tables the pipeline needs, at run
# time, into ${WORK}/metadata/. (No data/metadata files are stored in git.)
#   feature_ref.csv : Cell Ranger feature reference (TotalSeq-B hashtags, run MJZ019)
#   hash_map.csv    : hashtag -> mouse (JZ id) -> treatment arm
# Source of truth for any real run: 1-mouse_based/Tracing_Barcodes_list.xlsx
# =============================================================================
source config/config.sh
mkdir -p ${WORK}/metadata

# MJZ019 pooled 5 mice with hashtags B301/B302/B303/B304/B306 (B305 unused in this run).
cat > ${WORK}/metadata/feature_ref.csv <<'EOF'
id,name,read,pattern,sequence,feature_type
B301,B301,R2,5PNNNNNNNNNN(BC),ACCCACCAGTAAGAC,Antibody Capture
B302,B302,R2,5PNNNNNNNNNN(BC),GGTCGAGAGCATTCA,Antibody Capture
B303,B303,R2,5PNNNNNNNNNN(BC),CTTGCCGCATGTCAT,Antibody Capture
B304,B304,R2,5PNNNNNNNNNN(BC),AAAGCATTCTTCACG,Antibody Capture
B306,B306,R2,5PNNNNNNNNNN(BC),TATGCTGCCACGGTA,Antibody Capture
EOF

cat > ${WORK}/metadata/hash_map.csv <<'EOF'
hash_id,sample_id,model,treatment
B301,JZ201,Npp53,Castration + Enzalutamide
B302,JZ202,Npp53,Castration + Enzalutamide
B303,JZ203,Npp53,Castration + Enzalutamide
B304,JZ204,Npp53,Castration
B306,JZ205,Npp53,Castration
EOF
echo "wrote ${WORK}/metadata/{feature_ref,hash_map}.csv"
