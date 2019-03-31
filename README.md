<i>NB - this code / pipeline structure remains mostly untouched since it was first written by me in 2013 / 14 while working in the National Health Service (NHS) England - it would be designed differently were it written today. The main areas where it could be developed are:
  <ol><li>improved command line parameter parsing</li>
    <li>improved error handling, for example, checking output of each code section and making a decision to proceed or not</li></ol>

Placing the pipeline here on GitHub is primarily to archive it. Although it runs AOK in its current state, anyone re-using this code should make changes as you see fit.</i>

I update this pipeline as I re-use it myself in order to keep it maintained on a low level in line with new program versions that are released.

<b>Last update: Sunday, 31st March, 2019 @ 00:41 GMT</b>

# ClinicalGradeDNAseq - random read sampling to improve sensitivity
Automated next generation DNA sequencing analysis pipeline 'suited' for clinical tests, with >99.9% sensitivity to Sanger sequencing for <b><i>Single Nucleotide Variants</i></b> (SNVs) at read-depth>18 over target regions over interest.

This clinical-grade analysis pipeline, <i>ClinicalGradeDNAseq</i>, is a watered-down and modified version of the work by Blighe, Beauchamp, and colleagues at Sheffield Diagnostic Genetics Service, Sheffield Children's NHS Foundation Trust, Sheffield, UK, and their efforts to introduce a clinical-grade next generation sequencing (NGS) analysis pipeline fully validated against Sanger di-deoxy sequencing.

The pipeline is built using open source programs mixed with customised scripts. A wrapper script manages command line parameters and then executes the master analysis script. Control is then returned to the wrapper, where results files are transferred to a remote server via SSH/sFTP. A master and concise log is kept, with date- and time-stamps. Results directory structure is formed based on the run number and patient ID.

<b>[Key point]</b> The unique feature of the analysis pipeline that increases sensitivity to Sanger sequencing is in the variant calling step, where a final aligned BAM is split into 3 'sub-BAMs' of 75%, 50%, and 25% random reads. Variants are then called on all 4 BAMs, after which a consensus VCF is produced.

This pipeline has been modified from the original version and proceeds in an 8-step process:
<ol type="1">
  <li>Adaptor and read quality trimming - TrimGalore! (Krueger F), FastQC (Andrews S), cutadapt (Martin M, 2011)</li>
<li>Alignment - bwa mem (Li & Durbin, 2009)</li>
<li>Marking and removing PCR duplicates - Picard (Broad Institute of MIT and Harvard), SAMtools (Li et al., 2009)</li>
<li>Remove low mapping quality reads - SAMtools (Li et al., 2009)</li>
<li>QC - SAMtools (Li et al., 2009), custom scripts</li>
<li>Downsampling / random read sampling - Picard (Broad Institute of MIT and Harvard)</li>
<li>Variant calling - SAMtools/BCFtools (Li et al., 2009)</li>
<li>Annotation - Variant Effect predictor (McLaren et al., 2016)</li>
</ol>
      
<h1>Requirements</h1>
<ul>
  <li>Paired-end FASTQ files</li>
  <li>Chromosome-ordered and position-sorted BED file in hg19 / hg38 (BED files can be sorted with sort -k1,1V -k2,2n)</li>
  <li>Global installations of Cutadapt and unix2dos (included in dos2unix)</li>
  <li>Valid credentials for returning results files to remote server (username / password) via SSH/sFTP</li>
</ul>

<h1>Execution</h1>
<ol type="1">
Run the ‘PipelineWrapper’ wrapper script, which will check command-line parameters, execute the master script, and then return results files via SSH/sFTP to remote server. Use the following parameters:
<ol type="i"">
<li>FASTQ mate-pair 1 (absolute file path)</li>
<li>FASTQ mate-pair 2 (absolute file path)</li>
<li>Reference genome FASTA (absolute file path)</li>
<li>Run number (e.g. Plate6, Plate7, etc.) (alphanumeric)</li>
<li>Patient ID (alphanumeric)</li>
<li>BED file (absolute file path)</li>
<li>Minimum quality for bases at read ends, below which bases will be cut (integer)</li>
<li>Minimum allowed read length (integer)</li>
<li>Adaptor for trimming off read ends ('illumina' / 'nextera' / 'small_rna')</li>
<li>Minimum read depth for calling a variant (integer)</li>
<li>Minimum allowed mapping quality (integer)</li>
<li>Stringency for calling variants ('relaxed' / 'normal') (relaxed uses <code>--pval-threshold 1.0</code></li>
<li>Directory where results will be output (absolute file path)</li>
<li>User initials (alphanumeric)</li>
</ol>

<h1>Output</h1>
Results files are output locally to <i>[results root]/[run number]/[sample ID]/</i>, with a copy being also sent via SSH/sFTP to a remote server.
<ul>
  <li><i>*_AnalysisLog.txt</i> - analysis log (short)</li>
<li><i>Master.log</i> - analysis log (comprehensive)</li>
<li><i>*_R1_001.fastq.gz_trimming_report.txt</i> - details on base and read trimming for mate-pair 1</li>
<li><i>*_R1_001_val_1_fastqc.html</i> - FastQC report for mate-pair 1 (after trimming)</li>
<li><i>*_R2_001.fastq.gz_trimming_report.txt</i> - details on base and read trimming for mate-pair 2</li>
<li><i>*_R2_001_val_2_fastqc.html</i> - FastQC report for mate-pair 2 (after trimming)</li>
<li><i>*_Alignment.txt</i> - alignment metrics</li>
<li><i>*_ReadsOffTarget.txt</i> - number of reads falling outside regions specified in BED file</li>
<li><i>*_PCRDuplicates.txt</i> - details on identified PCR duplicates</li>
<li><i>*_CoverageTotal.bedgraph</i> - coverage for all mapped locations (contiguous bases at same read depth are merged into regions)</li>
<li><i>*_MeanCoverageBED.bedgraph</i> - mean read depth for each region specified in supplied BED file</li>
<li><i>*_PerBaseDepthBED.bedgraph</i> - per base read depth for each base in each region specified in supplied BED file</li>
<li><i>*_Aligned_Sorted_PCRDuped_FiltMAPQ.bam</i> - aligned BAM file with sorted reads, PCR duplicates removed, and reads below mapping quality threshold removed</li>
<li><i>*_Aligned_Sorted_PCRDuped_FiltMAPQ.bam.bai</i> - index for above BAM file</li>
<li><i>*_Final.vcf</i> - final VCF file</li>
<li><i>*_AnnotationVEP.html</i> - HTML report of variant annotation, with consequences for all known transcript isoforms</li>
<li><i>*_AnnotationVEP.tsv</i> - as above but in tab-separated values (TSV) format</li>
</ul>

<h1>Hard-coded sections of code</h1>
<ul>
  <li>PipelineWrapper.sh, line 120: <i>/home/ubuntu/pipeline/AnalysisMasterVersion1.sh "${Read1}" "${Read2}" ...</i> - absolute path filename for AnalysisMasterVersion1.sh</li>
  <li>PipelineWrapper.sh, line 138: <i>remoteDir="/remote/SAMBA/share/"</i> - Remote server directory to which results files will be transferred via SSH/sFTP</li>
  <li>PipelineWrapper.sh, line 150, 161: <i>sshpass -e sftp $username@XXX.XXX.XXX.XXX << !</i> - Remote server IP address or host name to which results files will be transferred via SSH/sFTP</li>
  <li>AnalysisMasterVersion1.sh, lines 25-35 - root directories (absolute paths) of required programs</li>
  <li>AnalysisMasterVersion1.sh, lines 286/304 - min-BQ set to 30 for BCFtools mpileup</li>
  <li>AnalysisMasterVersion1.sh, lines 374 - extra filters for filtering variants via vcfutils.pl (BCFtools)</li>
  <li>AnalysisMasterVersion1.sh, line 396/7 - 'species' and 'assembly' for VEP set to <i>homo_sapiens</i> and <i>GRCh38</i></li>
</ul>


<h1>References</h1>
<ul>
  <li>Andrews S, FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/, last accessed 28th August 2017.</li>
<li>Broad Institute of MIT and Harvard, Picard, http://broadinstitute.github.io/picard/, last accessed 28th August 2017</li>
<li>Krueger F, Trim Galore!, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/, last accessed 28th August 2017.</li>
<li>Li  H and Durbin R (2009), Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics 25(14): 1754–1760.</li>
<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup (2009), The Sequence Alignment/Map format and SAMtools, Bioinformatics 25(16):2078-9.</li>
<li>Martin M (2011), Cutadapt removes adapter sequences from high-throughput sequencing reads, EMBnet.journal 17(1): 10-12.</li>
<li>McLaren W, Gil L, Hunt S, Riat H, Ritchie G, Thormann A, Flicek P, Cunningham F (2016), The Ensembl Variant Effect Predictor, Genome Biology 17:122.</li>
</ul>
<h1>Credits</h1>
<ul>
  <li>Kevin Blighe (Sheffield Children's NHS Foundation Trust)</li>
  <li>Nick Beauchamp (Sheffield Children's NHS Foundation Trust)</li>
  <li>Darren Grafham (Sheffield Children's NHS Foundation Trust)</li>
  <li>Sheffield Diagnostics Genetics Service</li>
</ul>
