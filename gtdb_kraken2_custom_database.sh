cd /mnt/m3/xiaoning/rgcIntegr/kraken2_custom_database_R1000_R2048_bacteria/gtdb_to_taxdump

## remove 1st header line different from lineage2tax 
awk '(NR>1)' mag_gtdb_tax_formatted.tsv > mag_gtdb_tax_formatted_no_header.tsv

./gtdb_to_taxdump.py mag_gtdb_tax_formatted_no_header.tsv >log
 grep subspecies log > taxid_generate_mapping

## unified taxid for classification in header taxid derived from the  pseudo id generated from gtdb2tax 
mkdir seq_gtdb_v2
parallel  -j 36 --xapply 'sed "/>/ c\>{1}|kraken:taxid|{2}" mag_genome/{1}.fa >seq_gtdb_v2/{1}_pseudoTaxID.fa' ::: `cut -f 2 taxid_generate_mapping ` ::: `cut -f 1 taxid_generate_mapping`


## adding-library format or library-in-all-one format
mkdir ../bacteria_mag_gtdb_v2
parallel -j 36 'kraken2-build --add-to-library seq_gtdb_v2/{}_pseudoTaxID.fa --db ../bacteria_mag_gtdb_v2' ::: `cut -f 2 taxid_generate_mapping`



## prepare *.dmp file
mkdir -p ../bacteria_mag_gtdb_v2/taxonomy ../bacteria_mag_gtdb_v2/library
cp nodes.dmp names.dmp ../bacteria_mag_gtdb_v2/taxonomy/


## kraken2 build with CHT method(12min)
kraken2-build --build \
              --threads 48 \
              --db ../bacteria_mag_gtdb_v2

## pushd and popd
popd 
cp gtdb_to_taxdump/kraken2_custom_classify_gtdb_v2.sh . 
bash  kraken2_custom_classify_gtdb_v2.sh /mnt/m3/liufang/Rice_NRgene/01_host_remove_raw/kneaddata_output &



## say some batch integrate, how to update this custom database 
## for example integrate mag and isolate tsb

## step 0 gtdb tax file and seq preprocessing

ln -s /mnt/m3/liufang/Rice_NRgenome/yongxin_version/01_fasta/Rice_HQ isolate_tsb_genome
cp  /mnt/m3/liufang/Rice_NRgenome/yongxin_version/03_gtdb_tk/Rice/gtdbtk_output/Rice_isolate.bac120.summary.tsv isolate_tsb.bac120.tax.tsv
cut -f 1,2 isolate_tsb.bac120.tax.tsv  >isolate_tsb_tax_fmt.tsv

## step 1 concact the metadata without the header 
## metadata format: Genomeaccession\t d__B...;p__F...;...;s__B... at least 2 columns are mandetory

awk '(NR>1)' isolate_tsb_tax_fmt.tsv > isolate_tsb_tax_fmt_no_header.tsv

cat isolate_tsb_tax_fmt_no_header.tsv mag_gtdb_tax_formatted_no_header.tsv >integrate_metadata_tax.tsv

## step 2 generate ncbi *.dmp and pseudo id(order cannot be confirm each round)
./gtdb_to_taxdump.py integrate_metadata_tax.tsv >log_integrate

## extract mapping file mapped id to genomeaccession 
grep subspecies log_integrate > integrate_taxid_generate_mapping


 ## step 3 new library to build, may there be warning for part of file due to the whole mapping metadata
mkdir seq_gtdb_itegrt_v1
parallel  -j 36 \
          --xapply \
          'sed "/>/ c\>1{}|kraken:taxid|{2}" mag_genome/{1}.fa >seq_gtdb_itegrt_v1/{1}_pseudoTaxID.fa' \
          ::: `awk 'NR==FNR{a[$1]; next} ($1 in a){print }' <(cut -f 2 integrate_taxid_generate_mapping) <(cut -f 1 mag_gtdb_tax_formatted_no_header.tsv)` \
          ::: `awk 'NR==FNR{a[$1]; next} ($1 in a){print }' <(cut -f 2 integrate_taxid_generate_mapping) <(cut -f 1 mag_gtdb_tax_formatted_no_header.tsv)`


## warning .fa and .fna
## sed: can't read isolate_tsb_genome/R301.fna: No such file or directory
# sed: can't read isolate_tsb_genome/R285.fna: No such file or directory
# sed: can't read isolate_tsb_genome/R2165.fna:

parallel  -j 36 \
          --xapply \
          'sed "/>/ c\>{1}|kraken:taxid|{2}" isolate_tsb_genome/{1}.fna >seq_gtdb_itegrt_v1/{1}_pseudoTaxID.fa' \
          ::: `awk 'NR==FNR{a[$1]; next} ($1 in a){print }' <(cut -f 2 integrate_taxid_generate_mapping) <(cut -f 1 isolate_tsb_tax_fmt_no_header.tsv)` \
          ::: `awk 'NR==FNR{a[$1]; next} ($1 in a){print }' <(cut -f 2 integrate_taxid_generate_mapping) <(cut -f 1 isolate_tsb_tax_fmt_no_header.tsv)`


 
## step 4 adding-library format or library-in-all-one format
mkdir ../bacteria_gtdb_itegrt_v1
parallel -j 36 'kraken2-build --add-to-library seq_gtdb_itegrt_v1/{}_pseudoTaxID.fa --db ../bacteria_gtdb_itegrt_v1' ::: `cut -f 2 integrate_taxid_generate_mapping`


## step5 prepare *.dmp file 
mkdir -p ../bacteria_gtdb_itegrt_v1/taxonomy ../bacteria_gtdb_itegrt_v1/library
cp nodes.dmp names.dmp ../bacteria_gtdb_itegrt_v1/taxonomy/


## kraken2 build with CHT method(12min)
kraken2-build --build \
              --threads 48 \
              --db ../bacteria_gtdb_itegrt_v1


bash kraken2_custom_classify_gtdb_itegrt_v1.sh /mnt/m3/liufang/Rice_NRgene/01_host_remove_raw/kneaddata_output &

## output in kraken2/kraken2_reads  name format dereplicate 
bash manip.sh 

# phylum="d__ p__ c__ o__ f__ g__ s__"

# for tax in $phylum;
# do
#         echo $tax
#         sed -i "s/$tax$tax/$tax/g" kraken2/kraken2_reads/taxonomy_count_bacteria_gtdb_itegrt_v1.txt
# done