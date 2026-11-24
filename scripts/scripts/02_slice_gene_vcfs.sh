H="ftp.1000genomes.ebi.ac.uk"
D="vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"

# 1. APOE (chr19:44905796-44909393) -> 140 variants
bcftools view -S yri_samples.txt --force-samples -r chr19:44905796-44909393 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.recalibrated_variants.vcf.gz" -Oz -o ibadan_apoe.vcf.gz && bcftools index -f ibadan_apoe.vcf.gz

# 2. ABCA7 (chr19:1040101-1065571) -> 1,704 variants
bcftools view -S yri_samples.txt --force-samples -r chr19:1040101-1065571 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.recalibrated_variants.vcf.gz" -Oz -o ibadan_abca7.vcf.gz && bcftools index -f ibadan_abca7.vcf.gz

# 3. BIN1 (chr2:127048023-127107288) -> 2,734 variants
bcftools view -S yri_samples.txt --force-samples -r chr2:127048023-127107288 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.vcf.gz" -Oz -o ibadan_bin1.vcf.gz && bcftools index -f ibadan_bin1.vcf.gz

# 4. TREM2 (chr6:41158506-41163186) -> 218 variants
bcftools view -S yri_samples.txt --force-samples -r chr6:41158506-41163186 "http://${H}/${D}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.vcf.gz" -Oz -o ibadan_trem2.vcf.gz && bcftools index -f ibadan_trem2.vcf.gz
