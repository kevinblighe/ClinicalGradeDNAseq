
#!/bin/bash ;

echo -e "\n\n\n" ;
echo "###########################################################################" ;
echo "#Pre-liminary step:  Variable and log initialisation#" ;
echo "###########################################################################" ;

Read1="${1}" ;
Read2="${2}" ;
Ref_FASTA="${3}" ;
RunNumber="${4}" ;
PatientID="${5}" ;
BEDfile="${6}" ;
TrimmingQualityReadEnds="${7}" ;
TrimmingReadLengthMin="${8}" ;
TrimmingAdaptor="${9}" ;
FilterReadDepthCutoff="${10}" ;
FilterMappingQualityCutOff="${11}" ;
CallingStringency="${12}" ;
Results_root="${13}"
User="${14}" ;

#Root directories
Pipeline_root="/home/ubuntu/pipeline/" ;
Fastqc_root="/home/ubuntu/programs/FastQC/" ;
#Cutadapt_root="/home/ubuntu/programs/cutadapt-1.14/" ; #installed globally
TrimGalore_root="/home/ubuntu/programs/trim_galore_v0.4.4/" ;
BWA_root="/home/ubuntu/programs/bwa.kit/" ;
SAMtools_root="/home/ubuntu/programs/samtools-1.5/" ;
Picard_root="/home/ubuntu/programs/picard/" ;
BEDTools_root="/home/ubuntu/programs/bedtools2/bin/" ;
BCFtools_root="/home/ubuntu/programs/bcftools/" ;
VCFtools_root="/home/ubuntu/programs/vcftools_0.1.13/bin/" ;
VEP_root="/home/ubuntu/programs/ensembl-vep/" ;

#Begin log file
echo "Beginning analysis script on `date`, run by "${User}" with the following parameters:" > "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t1\t"${Read1}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t2\t"${Read2}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t3\t"${Ref_FASTA}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t4\t"${RunNumber}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t5\t"${PatientID}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t6\t"${BEDfile}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t7\t"${TrimmingQualityReadEnds}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t8\t"${TrimmingReadLengthMin}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t9\t"${TrimmingAdaptor}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t10\t"${FilterReadDepthCutoff}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t11\t"${FilterMappingQualityCutOff}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t12\t"${CallingStringency}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t13\t"${Results_root}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
echo -e "\t14\t"${User}"" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

echo "Done." ;



echo -e "\n\n" ;
echo "#####################################################" ;
echo "#Analysis step 1:  Adaptor and read quality trimming#" ;
echo "#####################################################" ;
echo "Beginning analysis step 1 (adaptor and read quality trimming) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

export PATH=""${Fastqc_root}":$PATH" ;

echo `"${TrimGalore_root}"/trim_galore --qual "${TrimmingQualityReadEnds}" --gzip --length "${TrimmingReadLengthMin}" --"${TrimmingAdaptor}" --paired "${Read1}" "${Read2}" --output_dir "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/ --fastqc` ;

read -r -a TrimmedFiles <<< `find "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/* -name "*.fq.gz"` ;

echo "Done." ;



echo -e "\n\n" ;
echo "#############################" ;
echo "#Analysis step 2:  Alignment#" ;
echo "#############################" ;
echo "Beginning analysis step 2 (alignment) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

#BWA mem
echo `"${BWA_root}"bwa mem "${Ref_FASTA}" "${TrimmedFiles[0]}" "${TrimmedFiles[1]}" > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.sam` ;

echo "Done." ;



echo -e "\n\n" ;
echo "###############################################" ;
echo "#Analysis step 3:  Marking PCR duplicate reads#" ;
echo "###############################################" ;
echo "Beginning analysis step 3 (marking and removing PCR duplicates) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

echo `"${SAMtools_root}"samtools view -bS "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.sam > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.bam` ;

echo `"${SAMtools_root}"samtools sort "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned.bam -o "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted.bam` ;

java -jar "${Picard_root}"picard.jar MarkDuplicates INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted.bam OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.bam ASSUME_SORTED=true METRICS_FILE="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.txt VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ;

echo `"${SAMtools_root}"samtools index "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.bam` ;

echo `"${SAMtools_root}"samtools view -b -F 0x400 "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDupes.bam > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.bam` ;

echo `"${SAMtools_root}"samtools index "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.bam` ;

echo "Done." ;



echo -e "\n\n" ;
echo "####################################################" ;
echo "#Analysis step 4:  Remove low mapping quality reads#" ;
echo "####################################################" ;
echo "Beginning analysis step 4 (remove low mapping quality reads) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

echo `"${SAMtools_root}"samtools view -bq "${FilterMappingQualityCutoff}" "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.bam > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam` ;

echo `"${SAMtools_root}"samtools index "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam` ;

echo "Done." ;



echo -e "\n\n" ;
echo "######################" ;
echo "#Analysis step 5:  QC#" ;
echo "######################" ;
echo "Beginning analysis step 5 (QC) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

echo `"${SAMtools_root}"samtools flagstat "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Alignment.txt` ;

echo `"${BEDTools_root}"bedtools intersect -v -bed -abam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam -b "${BEDfile}" | wc -l > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_ReadsOffTarget.txt` ;

#Output depth of coverage for all regions in the BAM file
#Sequential positions at the same read depth are merged into a single region
echo `"${BEDTools_root}"bedtools genomecov -ibam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam -bga -split > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_CoverageTotal.bedgraph` ;

#Output the per base read depth for each region in the BED file
echo `"${BEDTools_root}"bedtools coverage -a "${BEDfile}" -b "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam -d > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PerBaseDepthBED.bedgraph` ;

#Output the mean depth of coverage for each region in the BED file
echo `"${BEDTools_root}"bedtools coverage -a "${BEDfile}" -b "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam -mean > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_MeanCoverageBED.bedgraph` ;

#Get percent genome covered
zero=$("${BEDTools_root}"bedtools genomecov -ibam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam -g "${Ref_FASTA}" -bga | awk '$4==0 {bpCountZero+=($3-$2)} {print bpCountZero}' | tail -1)
nonzero=$("${BEDTools_root}"bedtools genomecov -ibam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam -g "${Ref_FASTA}" -bga | awk '$4>0 {bpCountNonZero+=($3-$2)} {print bpCountNonZero}' | tail -1)
percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
echo -e "Number of bases at 0 read-depth:\t""${zero}" > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt ;
echo -e "Number of bases at >0 read-depth:\t""${nonzero}" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt ;
echo -e "Percent reference genome covered:\t""${percent}" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt ;

echo "Done." ;



echo -e "\n\n" ;
echo "#######################################################" ;
echo "#Analysis step 6:  Downsampling / random read sampling#" ;
echo "#######################################################" ;
echo "Beginning analysis step 6 (downsampling / random read sampling) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

echo `java -jar "${Picard_root}"picard.jar DownsampleSam INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam RANDOM_SEED=50 PROBABILITY=0.75 VALIDATION_STRINGENCY=SILENT` ;

echo `"${SAMtools_root}"samtools index "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam` ;

echo `java -jar "${Picard_root}"picard.jar DownsampleSam INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam RANDOM_SEED=50 PROBABILITY=0.5 VALIDATION_STRINGENCY=SILENT` ;

echo `"${SAMtools_root}"samtools index "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam` ;

echo `java -jar "${Picard_root}"picard.jar DownsampleSam INPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam OUTPUT="${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam RANDOM_SEED=50 PROBABILITY=0.25 VALIDATION_STRINGENCY=SILENT` ;

echo `"${SAMtools_root}"samtools index "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam` ;

echo "Done." ;



echo -e "\n\n" ;
echo "###################################" ;
echo "#Analysis step 7:  Variant calling#" ;
echo "###################################" ;
echo "Beginning analysis step 7 (variant calling) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

if [ $CallingStringency == "relaxed" ] ;
then
	echo `"${SAMtools_root}"samtools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --output-tags DP,AD -f "${Ref_FASTA}" --BCF "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam | "${BCFtools_root}"bcftools call --consensus-caller --variants-only --pval-threshold 1.0 -Ob > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf` ;
elif [ $CallingStringency == "normal" ] ;
then
	echo `"${SAMtools_root}"samtools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --output-tags DP,AD -f "${Ref_FASTA}" --BCF "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_75pcReads.bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_50pcReads.bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_25pcReads.bam | "${BCFtools_root}"bcftools call --multiallelic-caller --variants-only -Ob > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf` ;
fi

#Re-header the VCF
echo -e ""${PatientID}"" > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.tsv ;
echo -e ""${PatientID}"_75pc" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.tsv ;
echo -e ""${PatientID}"_50pc" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.tsv ;
echo -e ""${PatientID}"_25pc" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.tsv ;
echo `"${BCFtools_root}"bcftools reheader --samples "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/samples.tsv "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf.tmp` ;

cp -f "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf.tmp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf ;

#1st pipe, /Programs/bcftools-1.3.1/bcftools norm -Ou -m -any
#	Left-align and normalize indels; split multiallelic sites into multiple rows
#	-Ou, output in uncompressed format
#	-m For multialleles, split any of SNPs or InDels
#
#2nd pipe, /Programs/bcftools-1.3.1/bcftools norm -Ou -f
#	As above but will check if REF alleles match the reference
#	-Ou, output in uncompressed format
#	-f, specify reference FASTA
echo `"${BCFtools_root}"bcftools norm -Ou -m-any "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bcf | "${BCFtools_root}"bcftools norm -Ov -f "${Ref_FASTA}" > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf` ;

#Filter the variants further
#Options:
#	-Q INT	minimum RMS mapping quality for SNPs [10]
#	-d INT	minimum read depth [2]
#	-D INT	maximum read depth [10000000]
#	-a INT	minimum number of alternate bases [2]
#	-w INT	SNP within INT bp around a gap to be filtered [3]
#	-W INT	window size for filtering adjacent gaps [10]
#	-1 FLOAT	min P-value for strand bias (given PV4) [0.0001]
#	-2 FLOAT	min P-value for baseQ bias [1e-100]
#	-3 FLOAT	min P-value for mapQ bias [0]
#	-4 FLOAT	min P-value for end distance bias [0.0001]
#	-e FLOAT	min P-value for HWE (plus F<0) [0.0001]
#	-p		print filtered variants
echo "Applying further filtering to called variants..." ;
echo "Variants filtered out:"
echo `"${BCFtools_root}"bcftools view "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.vcf | "${BCFtools_root}"misc/vcfutils.pl varFilter -d "${FilterReadDepthCutoff}" -w 1 -W 3 -a 1 -1 0.05 -2 0.05 -3 0.05 -4 0.05 -e 0.05 -p > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf` ;

#Sort the VCF
echo `"${VCFtools_root}"vcf-sort "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ_FiltExtra.vcf --chromosomal-order > "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Final.vcf` ;

echo "Done." ;



echo -e "\n\n" ;
echo "##############################" ;
echo "#Analysis step 8:  Annotation#" ;
echo "##############################" ;
echo "Beginning analysis step 8 (annotation) on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

echo `"${VEP_root}"vep --input_file "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Final.vcf --format vcf --species homo_sapiens --assembly GRCh38 --cache --refseq --fasta "${Ref_FASTA}" --check_ref --offline --bam "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam --output_file "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.tsv --tab --stats_file "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.html --variant_class --sift b --polyphen b --nearest symbol --gene_phenotype --regulatory --numbers --domains --vcf_info_field VEP --hgvs --hgvsg --symbol --tsl --canonical --af --af_1kg --af_esp --af_exac` ;



echo -e "\n\n\n" ;
echo "############################################################" ;
echo "#Post-analysis step:  Moving files from temporary directory#" ;
echo "############################################################" ;
echo "Beginning post-analysis tidy-up on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

#Convert text files to Windows format
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.txt ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Alignment.txt ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_ReadsOffTarget.txt ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_CoverageTotal.bedgraph ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PerBaseDepthBED.bedgraph ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_MeanCoverageBED.bedgraph ;
unix2dos "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt ;

#Copy analysis files
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped_FiltMAPQ.bam.bai "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Aligned_Sorted_PCRDuped.txt "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_PCRDuplicates.txt ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Alignment.txt "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_ReadsOffTarget.txt "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_CoverageTotal.bedgraph "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PerBaseDepthBED.bedgraph "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_MeanCoverageBED.bedgraph "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_PercentGenomeCovered.txt "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_Final.vcf "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.tsv "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/"${PatientID}"_AnnotationVEP.html "${Results_root}"/"${RunNumber}"/"${PatientID}" ;

read -r -a TrimmingReports <<< `find "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/* -name "*_trimming_report.txt"` ;
cp "${TrimmingReports[0]}" "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${TrimmingReports[1]}" "${Results_root}"/"${RunNumber}"/"${PatientID}" ;

read -r -a TrimmingHTMLReports <<< `find "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/* -name "*_fastqc.html"` ;
unix2dos "${TrimmingHTMLReports[0]}" ;
unix2dos "${TrimmingHTMLReports[1]}" ;
cp "${TrimmingHTMLReports[0]}" "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
cp "${TrimmingHTMLReports[1]}" "${Results_root}"/"${RunNumber}"/"${PatientID}" ;

rm -R "${Results_root}"/"${RunNumber}"/"${PatientID}"/tmp/ ;

chmod 755 "${Results_root}"/"${RunNumber}"/"${PatientID}" ;
chmod 755 "${Results_root}"/"${RunNumber}" ;

echo "Done." ;

echo "Analysis script finished on `date`" >> "${Results_root}"/"${RunNumber}"/"${PatientID}"/"${PatientID}"_AnalysisLog.txt ;

exit 0 ;
