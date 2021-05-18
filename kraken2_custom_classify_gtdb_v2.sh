#!/bin/bash
#########################################################################
# File Name: kraken2_classify.sh
# Author: Xiao-Ning Tank Zhang
# Email: xnzhang@genetics.ac.cn
# Created Time: Sun Feb  7 22:55:56 2021
#########################################################################
set -e
 
echo $0 $1 $2

# DBNAME=archaea
# DBNAME=bacteria
DBNAME=bacteria_mag_gtdb_v2
# DBNAME=bacteria_isolate_tsb
# DBNAME=bacteria_null_add2test
# DBNAME=bacteria_add2test
# BATCH=bacteria
BATCH=bacteria_mag_gtdb_v2
# BATCH=bacteria_null_dd2test_mag
# BATCH=bacteria_add2test_mag
# BATCH=bacteria_addR1000R2048

# INPUT_1=$1
# INPUT_2=$2
# THREADS=$2
# OUTPUT=$3

# kraken2 --db $DBNAME $INPUT_1 $INPUT_2 \
        # --threads 64 \
        # --output result_original \
        # --paired \
        # --report \
        # --use-mpa-style

# time parallel -j 6 \
#              kraken2 --db $DBNAME \
#             --paired $1/{1}*.fq.gz \
#             --threads 24 \
#             --use-names \
#             --use-mpa-style \
#             --report-zero-counts \
#             --report kraken2/kraken2_reads/{1}_report \
#             --output kraken2/kraken2_reads/{1}_output >> log/kraken2_reads/kraken2.log 2>&1 \
#             ::: `cat metadata.txt |cut -f 1`

# kraken2-build --add-to-library R1000_tax.fa --db $DBNAME
# kraken2-build --add-to-library R1001_tax.fa --db $DBNAME
# kraken2-build --build --db $DBNAME --threads 12

# time memusg parallel -j 6 \
#              kraken2 --db $DBNAME \
#             --paired $1/{1}*.fq.gz \
#             --threads 24 \
#             --use-names \
#             --use-mpa-style \
#             --report-zero-counts \
#             --report kraken2/kraken2_reads/{1}_report_addR1000R1001\
#             --output kraken2/kraken2_reads/{1}_output_addR1000R1001>> log/kraken2_reads/kraken2_addR1000R1001.log 2>&1 \
#             ::: `cat metadata.txt |cut -f 1`

# here is newest snippet
# time memusg parallel -j 6 \
#              kraken2 --db $DBNAME \
#             --paired $1/{1}*_kneaddata_paird_*\
#             --threads 24 \
#             --use-names \
#             --use-mpa-style \
#             --report-zero-counts \
#             --report kraken2/kraken2_reads/{1}_report_$BATCH\
#             --output kraken2/kraken2_reads/{1}_output_$BATCH>> log/kraken2_reads/kraken2_$BATCH.log 2>&1 \
#             ::: `cat metadata_id`


time memusg parallel -j 6 \
             kraken2 --db $DBNAME \
            --paired $1/{1}*kneaddata_paired_* \
            --threads 24 \
            --use-names \
            --use-mpa-style \
            --report-zero-counts \
            --report kraken2/kraken2_reads/{1}_report_$BATCH\
            --output kraken2/kraken2_reads/{1}_output_$BATCH>> log/kraken2_reads/kraken2_$BATCH.log 2>&1 \
            ::: `cat metadata_id`
## parallel kraken2 reads summary 

time memusg parallel -j 16 \
    'sort kraken2/kraken2_reads/{1}_report_{2} |cut -f 2 |sed "1 s/^/{1}\n/" \
    > kraken2/kraken2_reads/{1}_count_{2}' \
    ::: `cat metadata_id` ::: `echo $BATCH` 

kraken2_header=`cat metadata_id| head -n1 `
# kraken2_header=`cat bac_gtdb_taxid_test_final.tsv |cut -f 1| head -n1 `

sort kraken2/kraken2_reads/${kraken2_header}_report_$BATCH |cut -f 1 | sed "1 s/^/Taxonomy\n/" > kraken2/kraken2_reads/0header_count_$BATCH 
paste kraken2/kraken2_reads/*count_$BATCH> kraken2/kraken2_reads/taxonomy_count_$BATCH.txt



# time memusg parallel -j 16 \
#     'sort kraken2/kraken2_reads/{1}_report_addR1000R1001 |cut -f 2 |sed "1 s/^/{1}\n/" \
#     > kraken2/kraken2_reads/{1}_count_addR1000R1001' \
#     ::: `cat metadata.txt |cut -f 1`

# sort kraken2/kraken2_reads/${kraken2_header}_report_addR1000R1001 |cut -f 1 | sed "1 s/^/Taxonomy\n/" > kraken2/kraken2_reads/0header_count_addR1000R1001
# paste kraken2/kraken2_reads/*addR1000R1001 > kraken2/kraken2_reads/taxonomy_count_R1000R1001.txt






# parallel -j 36 'sed "/>/ c\>{}|kraken:taxid|{}" seq/{}/{}_contigs.fasta > seq/{}/{}_contigs_kraken_id.fa' ::: `cat seq/SampleID_tax_tag `

# parallel -j 36 'kraken2-build --add-to-library seq/{}/{}_contigs_kraken_id.fa --db archaea' ::: `cat seq/SampleID_tax_tag `
# cat seq/SampleID_tax_tag | sed "s/^R//g" |sed "s/[a-z]$//g" > seq/SampleID_tax_tag_only_num
# cut -f 1 bac_gtdb_tax_unified.tsv  | sed "s/^[A-z]//g" |sed "s/[A-z]$//g"  |paste - bac_gtdb_tax_unified.tsv >bac_gtdb_taxid.tsv
# awk '{print $1"\t|\tMOCKID_"$2"\t|\t\t|\t"$3"\t|"}' bac_gtdb_taxid.tsv >> archaea/taxonomy/names.dmp
# sed -i  '/^[A-z]/d' archaea/taxonomy/nodes.dmp
#  awk '{print $1"\t|\t1\t|\tspecies\t|\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|"}' bac_gtdb_taxid.tsv >> archaea/taxonomy/nodes.dmp
#  parallel -j 36 --xapply 'sed "/>/ c\>{2}|kraken:taxid|{1}" seq/{2}/{2}_contigs.fasta > seq/{2}/{2}_contigs_kraken_taxid_num.fa' ::: `cut -f 1 bac_gtdb_taxid.tsv ` ::: `cut -f 2 bac_gtdb_taxid.tsv`
# parallel -j 36 'kraken2-build --add-to-library seq/{}/{}_contigs_kraken_taxid_num.fa --db archaea' ::: `cut -f 2 bac_gtdb_taxid.tsv`
