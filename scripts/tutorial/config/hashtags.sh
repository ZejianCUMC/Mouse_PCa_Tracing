#!/bin/bash
# =============================================================================
# hashtags.sh — writes the two small metadata tables the pipeline needs, at run
# time, into ${WORK}/metadata/. (No data/metadata files are stored in git.)
#   feature_ref.csv : Cell Ranger feature reference (TotalSeq-B hashtags of this run)
#   hash_map.csv    : hashtag -> mouse (JZ id) -> model -> treatment arm
# The table is chosen by ${RUN} (set in config/config.sh).
# Source of truth for any real run: 1-mouse_based/Tracing_Barcodes_list.xlsx
# =============================================================================
source config/config.sh
mkdir -p ${WORK}/metadata

case "${RUN}" in
MJZ019*)
# MJZ019 pooled 5 Npp53 mice with hashtags B301/B302/B303/B304/B306 (B305 unused in this run).
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
;;
MJZ008*)
# MJZ008 (pilot run) pooled 3 intact mice with hashtags B301-B303.
# JZ136 has the Npp53 genotype; JZ137 / JZ138 are Np53Rb1 (Pten/p53/Rb1).
cat > ${WORK}/metadata/feature_ref.csv <<'EOF'
id,name,read,pattern,sequence,feature_type
B301,B301,R2,5PNNNNNNNNNN(BC),ACCCACCAGTAAGAC,Antibody Capture
B302,B302,R2,5PNNNNNNNNNN(BC),GGTCGAGAGCATTCA,Antibody Capture
B303,B303,R2,5PNNNNNNNNNN(BC),CTTGCCGCATGTCAT,Antibody Capture
EOF
cat > ${WORK}/metadata/hash_map.csv <<'EOF'
hash_id,sample_id,model,treatment
B301,JZ136,Npp53,Intact
B302,JZ137,Np53Rb1,Intact
B303,JZ138,Np53Rb1,Intact
EOF
;;
*)
echo "hashtags.sh: no hashtag table for RUN=${RUN}; add one from Tracing_Barcodes_list.xlsx" >&2; exit 1
;;
esac
echo "wrote ${WORK}/metadata/{feature_ref,hash_map}.csv for ${RUN}"
