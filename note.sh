# hyphen as positional parameter
# redirect <() must follow immediately after otherwise not wok
# paste to merge even and odd line
# paste to merge cut outputs

paste   <(cut -f 1 mag_gtdb_tax.tsv | sed -e '1s/ID/taxID/') <(cut -f 1 --complement mag_gtdb_tax.tsv | sed -e 's/\t/;/g' -e '1s/K.*/classification/g') >mag_gtdb_tax_formatted.tsv

# python '\t'.join(list1,list2)
awk '(NR>1){print NR+3000000  "\t" $s}' mag_gtdb_tax_formatted.tsv >mag_gtdb_tax_formatted_pseudoTaxID.tsv
./lineage2taxid.py  mag_gtdb_tax_formatted.tsv  ../bacteria/taxonomy/names.dmp ../bacteria/taxonomy/nodes.dmp
## generate gtdb_tax_final_formatted.tsv

# inner_join according to sampleID
join -1 2 -2 1    <(sort -k 2 mag_gtdb_tax_formatted_pseudoTaxID.tsv) <(sort -k 1 gtdb_tax_final_formatted.tsv) | sed 's/ /\t/g'> taxGTDB2NCBI.tsv
# join -1 1 -2 2   <(sort -k 1 gtdb_tax_final_formatted.tsv) <(sort -k 2 mag_gtdb_tax_formatted_pseudoTaxID.tsv) > taxGTDB2NCBI.tsv

mkdir $HOME/$nh/rgcIntegr/kraken2_custom_database_R1000_R2048_bacteria/gtdb_to_taxdump/seq
parallel  -j 36 --xapply 'sed "/>/ c\>{1}|kraken:taxid|{2}" mag_genome/{1}.fa >seq/{1}_pseudoTaxID.fa' ::: `cut -f 1 taxGTDB2NCBI.tsv ` ::: `cut -f 2 taxGTDB2NCBI.tsv `

awk '{print $2"\t|\t"$7"\t|\t"$6"\t|\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|"}' taxGTDB2NCBI.tsv >> ../bacteria_add2test/taxonomy/nodes.dmp
awk '{print $2"\t|\t"$3"\t|\t\t|\tscientific name\t|"}' taxGTDB2NCBI.tsv >> ../bacteria_add2test/taxonomy/names.dmp

parallel -j 36 'kraken2-build --add-to-library seq/{}_pseudoTaxID.fa --db ../bacteria_add2test' ::: `cut -f 1 taxGTDB2NCBI.tsv`
kraken2-build --build --threads 48 --db bacteria_add2test
