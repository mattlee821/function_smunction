#!/bin/bash

#SBATCH --job-name=join-test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# cp data
cd /user/home/ml16847/001_projects/adiposity_metabolites_endometrial_cancer/data/metabolites/phenospd/
cp /projects/MRC-IEU/research/projects/ieu3/p5/001/working/results/unadjusted/females/Ace*imputed* /user/home/ml16847/001_projects/adiposity_metabolites_endometrial_cancer/data/metabolites/phenospd/
ls Ace*imputed* > filelist

# make snp list
zcat Acetate-phase2-female-unadj_imputed.txt.gz | head
zcat Acetate-phase2-female-unadj_imputed.txt.gz | awk '{print $1, $2, $3, $5, $6}' | tail -n +2 | sort -k2,2 -k3,3 > temp
awk '{print $1, $2 ":" $3, $4, $5}' temp > snp_list
rm temp

# create data for each metabolite
export REMOVE=-phase2-female-unadj_imputed.txt.gz
for file in `cat filelist`; do
    # Step 1: Unzip 
    gzip -d -c ${file} > ${file}.unzipped # unzip GWAS
    
    # Step 2: Extract columns 1, 2, 3, 11, 12 and remvoe header
    cut -f1,2,3,11,12 ${file}.unzipped | tail -n +2 > ${file}.unzipped.extracted_columns

    # Sort the extracted columns file based on columns 2, and 3
    sort -k2,2 -k3,3 "${file}.unzipped.extracted_columns" > "${file}.unzipped.extracted_columns.sorted"
    
    # conbined column 2 and 3
    awk '{print $1, $2 ":" $3, $4, $5}' ${file}.unzipped.extracted_columns.sorted > ${file}.unzipped.extracted_columns.sorted.combined

    # Step 3: Left join using column 2 of both files
awk '
    NR == FNR     {f2[$2] = $3 OFS $4; next}
    $2 in f2      {print $0, f2[$2]}
' ${file}.unzipped.extracted_columns.sorted.combined snp_list > ${file}.unzipped.joined
 
    # step 4: add col name
    colNAME=$(echo "${file}" | sed -e "s/$REMOVE$//" | tr _ .)
    sed -i "1i SNP CHR:POS ALLELE1 ALLELE0 ${colNAME}_b ${colNAME}_se" "${file}.unzipped.joined"

    # Step 5: Extract columns 5 and 6 and create a new file
    cut -f5,6 "${file}.unzipped.joined" > "${file}.columns"
    
    # clean up
    rm ${file}.unzipped.joined ${file}.unzipped.extracted_columns.sorted.combined ${file}.unzipped.extracted_columns.sorted ${file}.unzipped.extracted_columns ${file}.unzipped 

done


# cut beta and se columns for each metabolite and make one combined data frame
ls *columns > filelist
# Loop through each file
for file in `cat filelist`; do
    # Extract columns 5 and 6 from the current file and combine them
    columns=$(awk '{print $5, $6}' "$file")
    if [ -z "$combined_columns" ]; then
        combined_columns="$columns"
    else
        combined_columns=$(paste <(echo "$combined_columns") <(echo "$columns"))
    fi
done
echo "$combined_columns" | column -t > combined_data
rm *columns

# make snp list for joining
echo -e "SNP ID ALLELE1 ALLELE0\n$(cat snp_list)" > snp_list_header
tr ' ' '\t' < snp_list_header > snp_list_header_formatted
cut -f1,3- snp_list_header_formatted > snp_list

# join snp list and combined data
paste snp_list combined_data > phenospd_data
rm snp_list_header_formatted snp_list_header snp_list filelist combined_data


# phenospd script
export SCRIPT=/user/home/ml16847/tools/PhenoSpD/script/my_phenospd.r
export DATA=/user/home/ml16847/001_projects/adiposity_metabolites_endometrial_cancer/data/metabolites/phenospd/
export OUT=/user/home/ml16847/001_projects/adiposity_metabolites_endometrial_cancer/data/metabolites/phenospd/all

cd /user/home/ml16847/tools/PhenoSpD/script/
module load lang/r/4.3.0-gcc
Rscript ${SCRIPT} --sumstats ${DATA}phenospd_data3 --out ${OUT}








