# ClinicalGradeDNAseq
Automated next generation DNA sequencing analysis pipeline suited for clinical tests

This clinical-grade analysis pipeline, ClinicalGradeDNAseq, is based on the work by Blighe, Beauchamp, and colleagues (Blighe et al., 2014) at Sheffield Diagnostic Genetics Service, Sheffield Children's National Health Service (NHS) Foundation Trust, Sheffield, UK, and their efforts to introduce a clinical-grade next generation sequencing (NGS) analysis pipeline fully validated against Sanger di-deoxy sequencing.
This pipeline has been modified from the original version and proceeds in a 9-step process:
<ol type="1">
  <li>Adaptor and read quality trimming - TrimGalore! (Krueger F), FastQC (Andrews S), cutadapt (Martin M, 2011)</li>
<li>Alignment - bwa mem (Li & Durbin, 2009)</li>
<li>Marking and removing PCR duplicates - Picard (Broad Institute of MIT and Harvard), SAMtools (Li et al., 2009)</li>
<li>Remove low mapping quality reads - SAMtools (Li et al., 2009)</li>
<li>QC - SAMtools (Li et al., 2009), custom scripts</li>
<li>Downsampling / random read sampling - Picard (Broad Institute of MIT and Harvard)</li>
<li>Variant calling - SAMtools/BCFtools (Li et al., 2009)</li>
<li>Annotation - Variant Effect predictor (McLaren et al., 2016)</li>
<li>Customising VCF (modifying indels; setting ID field to unique identifier; marking low quality variants) - custom scripts</li>
</ol>

<h1>Required input</h1>
<ul>
  <li>Paired-end FASTQ files</li>
<li>Chromosome ordered BED file in hg19 / hg38</li>
<li>The first column within all BED files must be in the format ‘chr1’, ‘chr22’, etc. BED files can be sorted with sort -k1,1V -k2,2n</li>
</ul>

<h1>Execution</h1>
<ol type="1">
Run the ‘PipelineKevin’ wrapper script, which will check command-line parameters, execute the master script, and then return results files via sFTP to remote server (password must be supplied). Use the following parameters:
<ol type="i"">
<li>FASTQ mate-pair 1 (absolute file path)</li>
<li>FASTQ mate-pair 2 (absolute file path)</li>
<li>Reference genome FASTA (absolute file path)</li>
<li>Plate number (e.g. Placa6, Placa7, etc.) (alphanumeric)</li>
<li>Patient ID (alphanumeric)</li>
<li>BED file (absolute file path)</li>
<li>PCR results CSV file (absolute file path)</li>
<li>PCR results reference patient ID (alphanumeric)</li>
<li>PCR results relative copy number gain cutoff (float)</li>
<li>PCR results relative copy number loss cutoff (float)</li>
<li>Minimum quality for bases at read ends, below which bases will be cut (integer)</li>
<li>Minimum allowed read length (integer)</li>
<li>Adaptor for trimming off read ends ('illumina' / 'nextera' / 'small_rna')</li>
<li>Minimum read depth for calling a variant (integer)</li>
<li>Minimum allowed mapping quality (integer)</li>
<li>Stringency for calling variants ('relaxed' / 'normal')</li>
<li>Directory where results will be output (absolute file path)</li>
<li>User initials (alphanumeric)</li>
             </ol>

<h1>Output</h1>
Results files are output locally to <i>[results root]/[run number]/[sample ID]/</i>, with a copy being also sent via sFTP to a remote server.
<ul>
  <li>*_AnalysisLog.txt - analysis log (short)</li>
<li>Master.log - analysis log (comprehensive)</li>
<li>*_R1_001.fastq.gz_trimming_report.txt - details on base and read trimming for mate-pair 1</li>
<li>*_R1_001_val_1_fastqc.html - FastQC report for mate-pair 1 (after trimming)</li>
<li>*_R2_001.fastq.gz_trimming_report.txt - details on base and read trimming for mate-pair 2</li>
<li>*_R2_001_val_2_fastqc.html - FastQC report for mate-pair 2 (after trimming)</li>
<li>*_Alignment.txt - alignment metrics</li>
<li>*_ReadsOffTarget.txt - number of reads falling outside regions specified in BED file</li>
<li>*_PCRDuplicates.txt - details on identified PCR duplicates</li>
<li>*_CoverageTotal.bedgraph - coverage for all mapped locations (contiguous bases at same read depth are merged into regions)</li>
<li>*_MeanCoverageBED.bedgraph - mean read depth for each region specified in supplied BED file</li>
<li>*_PerBaseDepthBED.bedgraph - per base read depth for each base in each region specified in supplied BED file</li>
<li>*_Aligned_Sorted_PCRDuped_FiltMAPQ.bam - aligned BAM file with sorted reads, PCR duplicates removed, and reads below mapping quality threshold removed</li>
<li>*_Aligned_Sorted_PCRDuped_FiltMAPQ.bam.bai - index for above BAM file</li>
<li>*_Final.vcf - final VCF file, which contains all variants that passed QC, and low quality ones marked as such</li>
<li>*_AnnotationVEP.html - HTML report of variant annotation, with consequences for all known transcript isoforms</li>
<li>*_AnnotationVEP.tsv - as above but in tab-separated values (TSV) format</li>
</ul>

<h1>Hard-coded sections of code</h1>


<h1>References</h1>
<ul>
  <li>Andrews S, FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/, last accessed 28th August 2017.</li>
<li>Blighe K, Beauchamp N, Allen KE, Nesbitt IM, Dawe J, Grafham D, Dalton A (2014), Next Generation Sequencing in the National Health Service England: A Pipeline that Completely Agrees with Sanger, Journal of Cancer Science and Therapy 6: 406-410.</li>
<li>Broad Institute of MIT and Harvard, Picard, http://broadinstitute.github.io/picard/, last accessed 28th August 2017</li>
<li>Krueger F, Trim Galore!, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/, last accessed 28th August 2017.</li>
<li>Li  H and Durbin R (2009), Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics 25(14): 1754–1760.</li>
<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup (2009), The Sequence Alignment/Map format and SAMtools, Bioinformatics 25(16):2078-9.</li>
<li>Martin M (2011), Cutadapt removes adapter sequences from high-throughput sequencing reads, EMBnet.journal 17(1): 10-12.</li>
<li>McLaren W, Gil L, Hunt S, Riat H, Ritchie G, Thormann A, Flicek P, Cunningham F (2016), The Ensembl Variant Effect Predictor, Genome Biology 17:122.</li>
</ul>
