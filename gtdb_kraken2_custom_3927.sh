# cd /mnt/m3/xiaoning/rgcIntegr/kraken2_custom_database_R1000_R2048_bacteria/gtdb_to_taxdump

# ## remove 1st header line different from lineage2tax 
# awk '(NR>1)' mag_gtdb_tax_formatted.tsv > mag_gtdb_tax_formatted_no_header.tsv

# ./gtdb_to_taxdump.py mag_gtdb_tax_formatted_no_header.tsv >log
#  grep subspecies log > taxid_generate_mapping

# ## unified taxid for classification in header taxid derived from the  pseudo id generated from gtdb2tax 
# mkdir seq_gtdb_v2
# parallel  -j 36 --xapply 'sed "/>/ c\>{1}|kraken:taxid|{2}" mag_genome/{1}.fa >seq_gtdb_v2/{1}_pseudoTaxID.fa' ::: `cut -f 2 taxid_generate_mapping ` ::: `cut -f 1 taxid_generate_mapping`


# ## adding-library format or library-in-all-one format
# mkdir ../bacteria_mag_gtdb_v2
# parallel -j 36 'kraken2-build --add-to-library seq_gtdb_v2/{}_pseudoTaxID.fa --db ../bacteria_mag_gtdb_v2' ::: `cut -f 2 taxid_generate_mapping`



# ## prepare *.dmp file
# mkdir -p ../bacteria_mag_gtdb_v2/taxonomy ../bacteria_mag_gtdb_v2/library
# cp nodes.dmp names.dmp ../bacteria_mag_gtdb_v2/taxonomy/


# ## kraken2 build with CHT method(12min)
# kraken2-build --build \
#               --threads 48 \
#               --db ../bacteria_mag_gtdb_v2

# ## pushd and popd
# popd 
# cp gtdb_to_taxdump/kraken2_custom_classify_gtdb_v2.sh . 
# bash  kraken2_custom_classify_gtdb_v2.sh /mnt/m3/liufang/Rice_NRgene/01_host_remove_raw/kneaddata_output &



## say some batch integrate, how to update this custom database 
## for example integrate mag and isolate tsb

## step 0 gtdb tax file and seq preprocessing

ln -s /mnt/m1/yongxin/rice/MAG/genome/IsolateMAG/temp/drep/0.95/dereplicated_genomes/ integrt_all_genome
cp  /mnt/m1/yongxin/rice/MAG/genome/IsolateMAG/result/gtdbtk/0.95/tax_ar.txt integrt_all_genome.tax.tsv

cut -f 1 integrt_all_genome.tax.tsv | paste -  <(cut --complement -f 1 integrt_all_genome.tax.tsv | sed 's/\t/;/g' ) >isolate_tsb_tax_fmt.tsv
# cut -f 1,2 integrt_all_genome.tax.tsv  >isolate_tsb_tax_fmt.tsv

## step 1 concact the metadata without the header 
## metadata format: Genomeaccession\t d__B...;p__F...;...;s__B... at least 2 columns are mandetory
# awk '(NR>1)' isolate_tsb_tax_fmt.tsv | cat - <(awk '(NR>1)' /mnt/m1/yongxin/rice/MAG/genome/all/temp/gtdbtk/g.ar122.summary.tsv | cut -f 1,2) > integrate_metadata_tax.tsv
awk 'NR>1' isolate_tsb_tax_fmt.tsv > integrate_metadata_tax.tsv
## check number of genomes
wc -l integrate_metadata_tax.tsv
## 1204 ok adding archeae



# cat isolate_tsb_tax_fmt_no_header.tsv mag_gtdb_tax_formatted_no_header.tsv >integrate_metadata_tax.tsv

## step 2 generate ncbi *.dmp and pseudo id(order cannot be confirm each round)
./gtdb_to_taxdump.py integrate_metadata_tax.tsv >log_integrate

## extract mapping file mapped id to genomeaccession 
grep subspecies log_integrate > integrate_taxid_generate_mapping


 ## step 3 new library to build, may there be warning for part of file due to the whole mapping metadata 
 ### check whether .fa or .fna
# mkdir seq_gtdb_itegrt_v1
# parallel  -j 36 \
#           --xapply \
#           'sed "/>/ c\>{1}|kraken:taxid|{2}" mag_genome/{1}.fa >seq_gtdb_itegrt_v1/{1}_pseudoTaxID.fa' \
#           ::: `awk 'NR==FNR{a[$1]; next} ($2 in a){print }' <(cut -f 1 mag_gtdb_tax_formatted_no_header.tsv) <(cat integrate_taxid_generate_mapping) |cut -f 2`\
#           ::: `awk 'NR==FNR{a[$1]; next} ($2 in a){print }' <(cut -f 1 mag_gtdb_tax_formatted_no_header.tsv) <(cat integrate_taxid_generate_mapping) |cut -f 1`
mkdir seq_gtdb_itegrt_v2
export LC_ALL='C'
parallel  -j 36 \
          --xapply \
          'sed "/>/ c\>{1}|kraken:taxid|{2}" integrt_all_genome/{1}.fna >seq_gtdb_itegrt_v2/{1}_pseudoTaxID.fa' \
          ::: `awk 'NR==FNR{a[$1]; next} ($2 in a){print }' <(cut -f 1 integrate_metadata_tax.tsv) <(cat integrate_taxid_generate_mapping) |cut -f 2`\
          ::: `awk 'NR==FNR{a[$1]; next} ($2 in a){print }' <(cut -f 1 integrate_metadata_tax.tsv) <(cat integrate_taxid_generate_mapping) |cut -f 1`


## warning .fa and .fna
## sed: can't read isolate_tsb_genome/R301.fna: No such file or directory
# sed: can't read isolate_tsb_genome/R285.fna: No such file or directory
# sed: can't read isolate_tsb_genome/R2165.fna:

# parallel  -j 36 \
#           --xapply \
#           'sed "/>/ c\>{1}|kraken:taxid|{2}" isolate_tsb_genome/{1}.fna >seq_gtdb_itegrt_v1/{1}_pseudoTaxID.fa' \
#           ::: `awk 'NR==FNR{a[$1]; next} ($2 in a){print }' <(cut -f 1 isolate_tsb_tax_fmt_no_header.tsv) <(cat integrate_taxid_generate_mapping) |cut -f 2`\
#           ::: `awk 'NR==FNR{a[$1]; next} ($2 in a){print }' <(cut -f 1 isolate_tsb_tax_fmt_no_header.tsv) <(cat integrate_taxid_generate_mapping) |cut -f 1`



 
## step 4 adding-library format or library-in-all-one format
mkdir ../bacteria_gtdb_itegrt_all_v1
parallel -j 36 'kraken2-build --add-to-library seq_gtdb_itegrt_v2/{}_pseudoTaxID.fa --db ../bacteria_gtdb_itegrt_all_v1' \
                ::: `cut -f 2 integrate_taxid_generate_mapping`


## step5 prepare *.dmp file 
mkdir -p ../bacteria_gtdb_itegrt_all_v1/taxonomy ../bacteria_gtdb_itegrt_all_v1/library
cp nodes.dmp names.dmp ../bacteria_gtdb_itegrt_all_v1/taxonomy/



## kraken2 build with CHT method(25min)
kraken2-build --build \
              --threads 48 \
              --db ../bacteria_gtdb_itegrt_all_v1


bash kraken2_custom_classify_gtdb_itegrt_v1.sh /mnt/m3/liufang/Rice_NRgene/01_host_remove_raw/kneaddata_output &

## output in kraken2/kraken2_reads  name format dereplicate 
bash manip.sh 

# phylum="d__ p__ c__ o__ f__ g__ s__"

# for tax in $phylum;
# do
#         echo $tax
#         sed -i "s/$tax$tax/$tax/g" kraken2/kraken2_reads/taxonomy_count_bacteria_gtdb_itegrt_v1.txt
# done


## tax summary and plot through log file 
#warning: log may append so each time to truncate the 212 line
# custom integr version
grep 'processed in' -A 1 log/kraken2_reads/kraken2_bacteria_gtdb_itegrt_v1.log  | grep -v '--' - | sed '$!N;s/\n/,/'|tail -n 212 > extr_cf_rt_bacteria_gtdb_itegrt_v1
# sed '$!N;s/\n/,/' extr_cf_rt_bacteria_gtdb_itegrt_v1

# custom mag only version
grep 'processed in' -A 1 log/kraken2_reads/kraken2_bacteria_mag_gtdb_v2.log   | grep -v '--' - | sed '$!N;s/\n/,/' |tail -n 212> extr_cf_rt_bacteria_mag_gtdb_v2
# sed '$!N;s/\n/,/' extr_cf_rt_bacteria_gtdb_itegrt_v1

# kraken2 standard bacteria version
grep 'processed in' -A 1 log/kraken2_reads/kraken2_bacteria_standard.log   | grep -v '--' - | sed '$!N;s/\n/,/' > extr_cf_rt_bacteria_standard



## kraken2 integrt all version
grep 'processed in' -A 1 log/kraken2_reads/kraken2_bacteria_gtdb_itegrt_all_v1.log  | grep -v '--' - | sed '$!N;s/\n/,/'|tail -n 212  > extr_cf_rt_bacteria_itegrt_all_v1

# inner join according the total reads number
join -1 1 -2 1  <(sort -k 1 extr_cf_rt_bacteria_standard ) <(sort -k 1 extr_cf_rt_bacteria_mag_gtdb_v2) > inner_join_extr_cf_std_mag
join -1 1 -2 1  <(sort -k 1 inner_join_extr_cf_std_mag ) <(sort -k 1 extr_cf_rt_bacteria_gtdb_itegrt_v1) > inner_join_extr_cf_std_mag_itergrt
join -1 1 -2 1  <(sort -k 1 inner_join_extr_cf_std_mag_itergrt) <(sort -k 1 extr_cf_rt_bacteria_itegrt_all_v1) > inner_join_extr_cf_std_mag_itergrt_all
cut -d ' ' -f 1,12,15,26,29,40,43,54,57 inner_join_extr_cf_std_mag_itergrt_all  | sed -e 's/ /\t/g' -e 's/(//g' -e 's/%)//g' >mag_3group_report.tsv

## mapping with linker of depth to sampleid
cut -f 1,6 final_kraken_mag_custom.tsv > sample_depth_mapping
join -1 2 -2 1 <(sort -k2 sample_depth_mapping) <(sort -k1 mag_3group_report.tsv) >  mag_custom_3group.tsv
awk -F ' ' '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10}' mag_custom_3group.tsv > a_final_kraken_mag_custom_3group.tsv

cp a_final_kraken_mag_custom_3group.tsv result/

cp /mnt/m3/liufang/Rice_NRgene/metadata/metadata_XN.csv .
sed  's/,/\t/g' metadata_XN.csv > metadata_xn.tsv
awk 'NR==FNR{a[$1]; next} ($1 in a){print }'  metadata_id metadata_xn.tsv >mapping_212.tsv
head -n 1 metadata_xn.tsv|cat - mapping_212.tsv >mappingfile_212.tsv 

## taxon distribution 






## diff check for ncbi and gtdb
# demo  
export LC_ALL="C"
grep '^C'  Y4211L_output_bacteria_mag_gtdb_v2   |cut -f 2|sort  > check_diff_custom_std/Y4211L_mag
grep '^C'  Y4211L_output_bacteria_cf_all_reads  |cut -f 2|sort > check_diff_custom_std/Y4211L_std

# comm  check_diff_custom_std/Y4211L_mag check_diff_custom_std/Y4211L_std
comm  -12 check_diff_custom_std/Y4211L_mag check_diff_custom_std/Y4211L_std | wc -l
## std unique 13; mag uniq 72; share 73

# diff check the overview across all the samples
time memusg parallel -j 6 \
            "grep '^C'  kraken2/kraken2_reads/{1}_output_bacteria_mag_gtdb_v2 | cut -f 2 | sort > kraken2/kraken2_reads/check_diff_custom_std/{1}_mag" \
            ::: `cat metadata_id`

time memusg parallel -j 6 \
            "grep '^C'  kraken2/kraken2_reads/{1}_output_bacteria_cf_all_reads | cut -f 2 | sort > kraken2/kraken2_reads/check_diff_custom_std/{1}_std" \
            ::: `cat metadata_id`

time memusg parallel -j 48 \
            "comm   kraken2/kraken2_reads/check_diff_custom_std/{1}_mag kraken2/kraken2_reads/check_diff_custom_std/{}_std | wc -l" \
            ::: `cat metadata_id`

# datamash transpose mapingfile_212 convenient for r plot
cut -d ' ' -f 1,2 result/a_final_kraken_mag_itegr_custom_3group_v1.tsv |sort -k1 |sed -i 's/ /\t/g'|datamash transpose >seq_depth



## output in kraken2/kraken2_reads  name format dereplicate such as p__p__
bash manip.sh 

# phylum="d__ p__ c__ o__ f__ g__ s__"

# change variable '-' into '.' for r processing meanwhile gtdb some bizzare names
sed  's/-/./g' taxonomy_count_bacteria_gtdb_itegrt_all_v1.txt >taxonomy_count_bacteria_gtdb_itegrt_all_v1_hypen2dot.txt
# sed  's/-/./g' taxonomy_count_bacteria_gtdb_itegrt_v1.txt >taxonomy_count_bacteria_gtdb_itegrt_v1_hypen2dot.txt

../../gtdb_to_taxdump/tax_rank_abundance_parse.py taxonomy_count_bacteria_gtdb_itegrt_all_v1_hypen2dot.txt -r 1 2 3 4 5 6 7 -o intgrt_all/
