params.threads = 4
params.outDir = "/workflow/output/"
params.beta = false
params.subsampling = true
params.max_depth = 80
params.min_depth = 30
params.max_bad_cov = 0.15
params.reads = "$baseDir/test/H37Rv_MB_R{1,2}.fastq.gz"
params.vcfplotter = "$baseDir/bin/vcfRplotter.R"
params.snpeff_base = "$baseDir/bin/snpEff"
params.centrifuge_ref_dir = "$baseDir/db/centrifuge_db"
params.centrifuge_ref  = "refv2"
params.centrifuge_cont = "$baseDir/db/centrifuge_db/refv2"
params.centrifuge_WGS = "$baseDir/db/centrifuge_db/refv2"
params.genome_db = "$baseDir/db/genomes"
params.resDB = "$baseDir/db/resistance_db"
params.hsp65 = "$baseDir/db/kma_hsp65/Tortoli_etal_hsp65"

// Parsing the input parameters
sampleName       = "$params.sampleName"
max_depth        = "$params.max_depth"
min_depth        = "$params.min_depth"
max_bad_cov      = "$params.max_bad_cov"
outDir           = "$params.outDir" + "/" + "$params.sampleName"
subsampling      = "$params.subsampling"
threads          = "$params.threads"
// databases
genome_db        = "$params.genome_db"
centrifuge_ref    = "$params.centrifuge_ref_dir" + "/" + "$params.centrifuge_ref"
conversion_table = "$params.centrifuge_ref_dir" + "/conversion_table.txt"
centrifuge_cont  = "$params.centrifuge_cont"
centrifuge_WGS   = "$params.centrifuge_WGS"
snpeff_base      = "$params.snpeff_base"
resDB            = "$params.resDB"
hsp65            = "$params.hsp65"
// beta
beta             = "$params.beta"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .set  { reads_ch1 }

// Tools paths and command prefixes
picard           = "java -jar $baseDir/bin/picard.jar"
snpit            = "/miniconda/envs/mycodentifier/bin/snpit" // https://github.com/philipwfowler/snpit
snpit_reporter   = "$baseDir/bin/MycoprofilerUtils/snpit_wrapper.py"
species_reporter = "$baseDir/bin/MycoprofilerUtils/species_report.py"
fastq_reporter   = "$baseDir/bin/MycoprofilerUtils/fastq_report.py"
merge_json       = "$baseDir/bin/MycoprofilerUtils/merge_json.py"
reporter         = "$baseDir/bin/MycoprofilerUtils/final_report.py"
KMA_reporter     = "$baseDir/bin/MycoprofilerUtils/KMA_report.py"
add_stats        = "$baseDir/bin/MycoprofilerUtils/add_stats.py"
check_depth      = "$baseDir/bin/MycoprofilerUtils/check_zero_depth.py"
convert_tbdb     = '$baseDir/bin/convert_tbdb.py'
consensus_vcf    = "$baseDir/bin/indel_consensus.py"
generate_db      = "$baseDir/bin/gen_snpeff_db.py"
snpsift          = "java -jar $baseDir/bin/SnpSift/SnpSift.jar"
resis_reporter   = "$baseDir/bin/resistance_reporter.py"
VCFplotter       = "$params.vcfplotter"
error_strategy   = { sleep(Math.pow(2, task.attempt) * 500 as long); return 'retry' }
max_retries      = 1

log.info """
NEXTFLOW MyCodentifier V1.0
================================
sampleName     : $params.sampleName
reads      : $params.reads
outDir     : $params.outDir
codeBase   : $baseDir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
threads    : $params.threads
subsampling: $params.subsampling
max_depth  : $params.max_depth
min_depth  : $params.min_depth
max_bad_cov: $params.max_bad_cov
beta       : $params.beta
~~~~~~~~~~~Databases~~~~~~~~~~~
genome       : $genome_db
references   : $centrifuge_ref
TB/NTM       : $centrifuge_WGS
contaminents : $centrifuge_cont
resDB        : $resDB
hsp65        : $hsp65
================================
"""

/*********************************************************************************
 * PART 1: Read cleaning, organism identification and fetching of reference genome
 *********************************************************************************/

// Clean reads (adapter and read length filter)
process '1A_clean_reads' {
 tag '1A'
 errorStrategy error_strategy
 maxRetries max_retries
 conda 'bioconda::fastp=0.20.0 bioconda::pyfastx=0.6.12 conda-forge::simplejson=3.17.0'
 publishDir outDir + '/json', mode: 'copy', pattern: "*.fastp_json"
 input:
   set pairID, file(reads) from reads_ch1
 output:
   set file("${reads[0].baseName}_fastp.fastq.gz"), file("${reads[1].baseName}_fastp.fastq.gz") into fastp_1B, fastp_2A, fastp_1E, fastp_1G, fastp_1I
   file "${sampleName}.fastp_json" into fastp_json
   file "${sampleName}_*.json" into meta_json_fastq
   file ".command.*"
 script:
   """
   $fastq_reporter -i ${reads[0]} --sample ${sampleName}
   
   fastp -i ${reads[0]} -I ${reads[1]} -o ${reads[0].baseName}_fastp.fastq.gz -O ${reads[1].baseName}_fastp.fastq.gz \
   --length_required 50 --json ${sampleName}.fastp_json --html ${sampleName}.fastp.html --thread ${threads}

   # merge json results to meta json
   $merge_json --json1 ${sampleName}_*.json --json2 ${sampleName}.fastp_json --key fastp
   """
}

// Process 1B: Identify mycobacterium organisms (mainly for metagenomic data)
process '1B_ID_strain' {
  tag '1B'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::centrifuge=1.0.4_beta'
  publishDir outDir + '/centrifuge', mode: 'copy'
  input:
  file reads from fastp_1B
  output:
    file "${sampleName}_report.file" into report_file
    file ".command.*"
  script:
    """
    if [ `ls -1 ${centrifuge_ref}* 2>/dev/null | wc -l ` -gt 0 ];
    then
      centrifuge -p ${threads} -x ${centrifuge_ref} --quiet -q -1 ${reads[0]} -2 ${reads[1]} \
      --report-file ${sampleName}_report.file > ${sampleName}cent_stdout.out
    else
      echo "Centrifuge db: ${centrifuge_ref} does not exist."
      exit 1
    fi
    """
}

// Process 1C: Read the report from ID_strain step and store the data in a JSON with custom module
process '1C_read_ID_strain' {
  tag '1C'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'conda-forge::jq=1.6'
  publishDir outDir + '/centrifuge', mode: 'copy'
  input:
    file report from report_file
    file meta from meta_json_fastq
  output:
    file meta into meta_json_contaminant
    stdout() into (acc_1D, acc_1G, acc_1H, acc_2B, acc_3A, acc_4E, acc_5A, acc_5B, acc_5C, acc_6A)
    file("${sampleName}_spreporter.out")
    file ".command.*"
  script:
    """
    $species_reporter --JSON ${meta} -i $report --sample ${sampleName} --conv_table ${conversion_table} > ${sampleName}_spreporter.out
    cat ${sampleName}_*.json | jq .centrifuge.accession | xargs echo -n
    """
}

// Process 1D: Copy the fasta, genbank (GFF) files from db to be used by processes
process '1D_copy_genome_files' {
    errorStrategy error_strategy
    maxRetries max_retries
    conda 'conda-forge::requests=2.23.0 conda-forge::parallel=20200322'
    tag '1D'
    input:
     val acc from acc_1D
    output:
        file "${acc}_filt_genomic.fa" into (genome_2A, genome_2B, genome_2D, genome_3A, genome_3E, genome_5A, genome_3D, genome_4A, genome_4B, genome_4D)
        file "${acc}_filt_genomic.fa.fai" into (genome_idx_2B, genome_idx_2D, genome_index3D, genome_idx_4A, genome_idx_4B, genome_idx_4D)
        file "${acc}_genomic.gff" into (gff_5A, gff_5B)
        file "${acc}_genomic.gbff" into gbk_5A
        file ".command.*"
    script:
    """
    genome_dir="${genome_db}/${acc}"
    
    if [ -d \${genome_dir} ]
    then
      cp \${genome_dir}/${acc}_filt_genomic.fa ${acc}_filt_genomic.fa
      cp \${genome_dir}/${acc}_filt_genomic.fa.fai ${acc}_filt_genomic.fa.fai
      cp \${genome_dir}/${acc}_genomic.gff ${acc}_genomic.gff
      cp \${genome_dir}/${acc}_genomic.gbff ${acc}_genomic.gbff
    else
      echo "Please build the genome database first."
      exit 1
    fi
    """
}

// Process 1E: Measure probable contaminants (map metagenomic WGS data to a database containing all organism except mycobacteria)
process '1E_contaminants' {
  tag '1E'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::centrifuge=1.0.4_beta'
  publishDir outDir + '/contaminants', mode: 'copy'
    input:
        file reads from fastp_1E
    output:
        file "${sampleName}_report.file" into report_file_1F
        file ".command.*"
    script:
    """
    if [ `ls -1 ${centrifuge_cont}* 2>/dev/null | wc -l ` -gt 0 ];
    then
      centrifuge -p ${threads} -x ${centrifuge_cont} --quiet -q -1 ${reads[0]} -2 ${reads[1]} \
      --report-file ${sampleName}_report.file --exclude-taxids 1763 > ${sampleName}cent_stdout.out
    else
      echo "Centrifuge db: ${centrifuge_cont} does not exist."
      exit 1
    fi
    """
}

// Process 1F: Measure probable contaminants (map metagenomic WGS data to a database containing all organism except mycobacteria)
process '1F_read_contaminants' {
  tag '1F'
  publishDir outDir + '/contaminants', mode: 'copy'
  input:
    file report from report_file_1F
    file meta from meta_json_contaminant
  output:
    file "${meta.baseName}.json" into meta_json_1H
    file ".command.*"
  script:
    """
    $species_reporter --JSON ${meta} -i $report --sample ${sampleName} --omit_accession \
    --name contaminants > ${sampleName}_spreporter.out
    """
}

// Process 1G: hsp65 database to identify NTM/TB
process '1G_hsp65_identification' {
  tag '1G'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::kma=1.3.2'
  publishDir outDir + '/hsp65', mode: 'copy'
  input:
    file reads from fastp_1G
    val acc from acc_1G
  output:
    file "hsp65.res" into hsp65res_1H
    file "hsp65.aln"
    file ".command.*"
  script:
    """
    kma -t_db ${hsp65} -ipe ${reads[0]} ${reads[1]} -o hsp65
    """
}

// Process 1H: read hsp65 results
process '1H_read_hsp65' {
  tag '1H'
  publishDir outDir + '/hsp65', mode: 'copy'
  input:
    file meta from meta_json_1H
    file results from hsp65res_1H
    val acc from acc_1H
  output:
    file "${meta.baseName}.json" into meta_json_1J
    file ".command.*"
  script:
    """
    $KMA_reporter -i ${results} --database hsp65 --sample ${sampleName} --JSON ${meta}
    """
}

// Process 1I: WGS identify with centrifuge NTM/TB
process '1I_WGS_typing' {
  tag '1I'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::centrifuge=1.0.4_beta'
  publishDir outDir + '/WGS_typing', mode: 'copy'
    input:
        file reads from fastp_1I
    output:
        file "${sampleName}_report.file" into report_file_1J
        file ".command.*"
    script:
    """
    if [ `ls -1 ${centrifuge_WGS}* 2>/dev/null | wc -l ` -gt 0 ];
    then
      centrifuge -p ${threads} -x ${centrifuge_WGS} --quiet -q -1 ${reads[0]} -2 ${reads[1]} \
      --report-file ${sampleName}_report.file --exclude-taxids 1763 > ${sampleName}cent_stdout.out
    else
      echo "Centrifuge db: ${centrifuge_WGS} does not exist."
      exit 1
    fi
    """
}

// Process 1J: Retrieve WGS Typing results
process '1J_read_WGS_typing' {
  tag '1J'
  publishDir outDir + '/WGS_typing', mode: 'copy'
  input:
    file report from report_file_1J
    file meta from meta_json_1J
  output:
    file "${meta.baseName}.json" into meta_json_2B
    file ".command.*"
  script:
    """
    $species_reporter --JSON ${meta} -i $report --sample ${sampleName} --omit_accession \
    --name WGS_typing > ${sampleName}_spreporter.out
    """
}


/*************************************************************
 * PART 2: Read mapping, BAM QC and some BAM post-processing
 *************************************************************/

// Process 2A: Map the reads to the reference genome
process '2A_map_paired_reads' {
  tag '2A'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.9'
  input:
     file genome from genome_2A
     set file(trimmed1), file(trimmed2) from fastp_2A
  output:
     file("${sampleName}.bam") into (bam_2B, bam_2C)
     file("${sampleName}.bam.bai") into (bam_idx_2B, bam_idx_2C)
     file ".command.*"
  script:
    """
     bwa index $genome
     bwa mem -M -t ${threads} $genome ${trimmed1} ${trimmed2} | samtools sort -@2 -o ${sampleName}.bam -
     samtools index ${sampleName}.bam > ${sampleName}.bam.bai
    """
}

// Process 2B: Do BAM QC, report QC in JSON
process '2B_BAM_qc' {
  tag '2B'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::goleft=0.2.0 conda-forge::jq=1.6'
  publishDir outDir + '/bam_qc', mode: 'copy'
  input:
     file bam from bam_2B
     file idx from bam_idx_2B
     val acc from acc_2B
     file meta from meta_json_2B
     file genome from genome_2B
     file genome_idx from genome_idx_2B
  output:
      file "${meta.baseName}_qc.json" into meta_json_resistance, meta_json_snpit, meta_json_2C, meta_json_nocov
      file "BAM_QC.txt"
      file "avg_dp.txt"
      file ".command.*"
  script:
   """
   # boolean value QC passed
   passed=1

   # compute coverage with 16kb resolution to do quick bamQC
   # override exclude pattern, default excludes: ^NC
   goleft indexcov --excludepatt ^no_excl\$ --directory data/coverage $bam

   # get the percentage of high/low coverage from output
   low=\$(zcat data/coverage/coverage-indexcov.bed.gz | awk '{ if (\$4 < 0.15) print ; }' | wc -l)
   high=\$(zcat data/coverage/coverage-indexcov.bed.gz | awk '{ if (\$4 < 0.85) print ; }' | wc -l)
   all=\$(zcat data/coverage/coverage-indexcov.bed.gz | wc -l)

    # if the precentage of 16kb chunks with low coverage reads is higher than a percentage (default: 15%), terminate pipeline
   low_cov_frac=\$(echo "scale=2; \$low / \$all" | bc)
   high_cov_frac=\$(echo "scale=2; \$high / \$all" | bc)
   if ((  \$( echo "\$low_cov_frac > $max_bad_cov" | bc -l ) ));
   then
      msg=("BAM coverage is too low! (percentage of low coverage reads: \${low_cov_frac} >  $max_bad_cov)")
      echo \$msg > BAM_QC.txt
      passed=0;
   else
     msg=("BAM QC PASSED! (percentage of low coverage reads: \${low_cov_frac} <  $max_bad_cov)")
     echo \$msg > BAM_QC.txt;
   fi

   # calculate average depth of bam file, in order to determine the percentage to throw away
   goleft depth $bam --reference $genome --prefix $sampleName --mincov 0 --windowsize 10

   avg=\$(cat "${sampleName}.depth.bed" | awk '{ sum += \$4 } END { if (NR > 0) print sum / NR }')
   # stop pipeline if depth is too low (minimum depth default=30)
   if (( \$( echo "\${avg} < $min_depth" |bc -l ) ));
   then
     msg=("Average depth too low! ( \${avg} < cutoff: $min_depth )")
     echo \$msg
     echo \$msg >> BAM_QC.txt
     passed=0;
   else
     msg=("Average depth passed threshold ( \$avg > cutoff: $min_depth )")
     echo \$msg
     echo \$msg >> BAM_QC.txt;
   fi;
   # write average to a text file for subsampling
   echo \$avg > avg_dp.txt

   # add QC info to meta json file
   qc='{"avg_dp": '"\$avg"' ,"low_cov_frac": '"\$low_cov_frac"', "dp_thres": '"$min_depth"' , "cov_thres": '"$max_bad_cov"' , "qc_passed": '"\$passed"' }'

   echo \$qc
   echo \$(jq --argjson addqc "\$qc" '. + {BAM_QC : \$addqc }' ${meta}) > ${meta.baseName}_qc.json
   echo  ${meta.baseName}_qc.json
   [ -s ${meta.baseName}_qc.json ] || echo "file is empty" || exit 1;
 """
}

// Process 2C: filter the unmapped reads from the BAM file and subsample the BAM (reduce the size for faster processing time) if needed
process '2C_subsample_BAM' {
  tag '2C'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::samtools=1.9 bioconda::sambamba=0.6.6 conda-forge::jq=1.6'
  input:
    file bam from bam_2C
    file bam_idx from bam_idx_2C
    file meta from meta_json_2C
  output:
      file "${bam.baseName}*_2C.bam" into bam_2D
      file "${bam.baseName}*_2C.bam.bai" into bam_idx_2D
      file ".command.*"
  script:
  if(subsampling == false )
    """
      # throw away unmapped reads, make new bam with corresponding index
      samtools view -b -F 4 -f 3 -q 60 $bam > ${bam.baseName}_mapped.bam
      samtools index ${bam.baseName}_mapped.bam > ${bam.baseName}_mapped.bam.bai

      mv ${bam.baseName}_mapped.bam ${bam.baseName}_mapped_2C.bam
      mv ${bam.baseName}_mapped.bam.bai ${bam.baseName}_mapped_2C.bam.bai
    """
  else if(subsampling == true)
    """
      # get avg depth from file
      avg=\$(jq '.BAM_QC.avg_dp' ${meta} | xargs echo -n)
      # calculate the percentage of reads to throw away when we want to keep an average of depth 60 (default)
      perc="\$(echo "scale=2; ${max_depth} / \${avg}" | bc)"
      echo \$avg
      echo \$perc
      # throw away unmapped reads, make new bam with corresponding index
      samtools view -b -F 4 -f 3 -q 60 $bam > ${bam.baseName}_mapped.bam
      samtools index ${bam.baseName}_mapped.bam > ${bam.baseName}_mapped.bam.bai

      # if the percentage is higher than 1, sambamba will not modify the file in the subsampling process
      echo "subsampling BAM, keeping \${perc} percentage of the data for an avarage read depth of: ${max_depth}"
      # subsample bam file, throw away a percentage of the reads to achieve an average read depth of user specified value, default=80
      sambamba view -h -t 4 -s \$perc -f bam --subsampling-seed=6 ${bam.baseName}_mapped.bam -o ${bam.baseName}_subsampled.bam
      samtools index ${bam.baseName}_subsampled.bam > ${bam.baseName}_subsampled.bam.bai

      mv ${bam.baseName}_subsampled.bam ${bam.baseName}_subsampled_2C.bam
      mv ${bam.baseName}_subsampled.bam.bai ${bam.baseName}_subsampled_2C.bam.bai
  """
}

// Process 2D: Filter bad mapped reads
process '2D_Filter_mapped_reads' {
  tag '2D'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'bioconda::samtools=1.9 bioconda::samclip=0.2 bioconda::bamutil=1.0.14'
  input:
     file genome from genome_2D
     file genome_idx from genome_idx_2D
     file bam from bam_2D
     file bam_idx from bam_idx_2D
  output:
     file "${sampleName}.final.bam" into bam_3B
     file "${sampleName}.final.bam.bai" into bam_idx_3B
     file ".command.*"
  script:
    """

    # perform bamUtils to remove and trim reads with many mismatches
    bam filter --in ${bam} --refFile ${genome} --out filter.bam
    samtools view -b -F 4 -f 3 -q 60 filter.bam > fixed.filter.bam

    $picard FixMateInformation TMP_DIR=${outDir}/tmp ASSUME_SORTED=false I=fixed.filter.bam O=filter.fixed.nsorted.bam

    samtools sort filter.fixed.nsorted.bam > filter.sorted.fixed.nsorted.bam
    samtools index filter.sorted.fixed.nsorted.bam

    # remove all soft and hardclipped reads
    samtools view -h filter.sorted.fixed.nsorted.bam | samclip --ref ${genome} > final.sam
    samtools view -Sb final.sam | samtools sort - > ${sampleName}.final.bam
    samtools index ${sampleName}.final.bam
    """
}

/**********************************************************************************************
 * PART 3: Do the preprocessing steps for variant calling with GATK and execute HaplotypeCaller
 **********************************************************************************************/

//Process 3A: Create a sequence dictionary for GATK
process '3A_create_SeqDictionary' {
  tag '3A'
  publishDir outDir + '/' + "${acc}", mode: 'copy'
  input:
      val acc from acc_3A
      file genome from genome_3A
  output:
      // genome_dict formatted like SAM header, describing the contents of fasta file
      file "${genome.baseName}.dict" into (genome_dict3D, fasta_dict_4B, fasta_dict_4A, fasta_dict_4D)
      file ".command.*"
  script:
    """
    $picard CreateSequenceDictionary R=${genome} O=${genome.baseName}.dict
    """
}

// Process 3B: Mark read groups (.._rg.bam)
process '3B_Read_Group_GATK' {
  tag '3B'
  publishDir outDir + '/bamfiles'
  input:
    file bam from bam_3B
    file bam_idx from bam_idx_3B
  output:
    file "${bam.baseName}_rg.bam" into bam_3C
    file "${bam.baseName}_rg.bam.bai" into bam_idx_3C
    file ".command.*"
  script:
    """
    # use parallel for this step if it proves to be resource intensive?
    $picard AddOrReplaceReadGroups I=${bam} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${bam.baseName} O=${bam.baseName}_rg.bam
    samtools index ${bam.baseName}_rg.bam > ${bam.baseName}_rg.bam.bai
    """
}

// Process 3C: Mark read groups (.._rg.bam)
process '3C_MarkDuplicates_GATK' {
  tag '3C'
  publishDir outDir + '/bamfiles'
  input:
    file bam from bam_3C
    file bam_idx from bam_idx_3C
  output:
    file "${bam.baseName}_MD.bam" into (bam_3D, bam_3E, bam_5C)
    file "${bam.baseName}_MD.bam.bai" into (bam_idx_3D, bam_idx_3E)
    file ".command.*"
  script:
    """
    # use parallel for this step if it proves to be resource intensive?
    $picard MarkDuplicates REMOVE_DUPLICATES=true I=${bam} O=${bam.baseName}_MD.bam M=marked_dup_metrics.txt
    samtools index ${bam.baseName}_MD.bam > ${bam.baseName}_MD.bam.bai
    """
}

// Process 3D: run HaplotypeCaller to call the variants
process '3D_runHaplotypeCaller' {
  tag '3D'
  publishDir outDir + '/GATK', mode: 'copy'
  input:
    file genome from genome_3D
    file fast_idx from genome_index3D
    file dict from genome_dict3D
    file bam from bam_3D
    file bam_idx from bam_idx_3D
  output:
    file "raw.snps.indels_${bam.baseName}.vcf" into GATK_4A, GATK_4B
    file ".command.*"
  script:
    // stand_call_conf: PHRED score cutoff
    // nct, nt: number of cores, number of threads per core (not available for this caller)
    // what is right value for ploidy??
    // -nct ${threads}
    // --bamOutput gatk_${bam.baseName}.bam
    // needed params for GVCF: -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000
    """
    ##GATK still add -stand_call_conf 30 to GATK4
    ##gatk --java-options "-Xmx4G" HaplotypeCaller --sample-ploidy 2 -R ${genome} -I ${bam} -O raw.snps.indels_${bam.baseName}.vcf

    gatk --java-options "-Xmx8G" HaplotypeCaller --sample-ploidy 2 -R ${genome} -I ${bam} -O raw.snps.indels_${bam.baseName}.vcf \
    --max-assembly-region-size 600 --standard-min-confidence-threshold-for-calling 30.0 --min-assembly-region-size 300
    """
}

// Process 3E: run Delly to call large structural variants
process '3E_Detect_large_SV' {
    errorStrategy error_strategy
    maxRetries max_retries
    conda 'bioconda::delly=0.7.6 bioconda::samtools=1.9'
    tag '3E'
    publishDir outDir + '/delly', mode: 'copy'
    input:
        file genome from genome_3E
        file bam from bam_3E
        file idx from bam_idx_3E
    output:
        file "${sampleName}_delly.vcf.gz" into delly_vcf
        file ".command.*"
    script:
    """
    delly call -t DEL -g ${genome} ${bam} -o ${sampleName}_delly.bcf
    # bcftools query -f "%CHROM\t%POS\t[%END\t%GT\t%DR\t%DV\t%RR\t%RV]\\n" data/bamfiles/ERR1664619.delly.bcf
    # compress and index vcf file

    # check if file exists (if not no deletions are found copy mock/empty file to continue)
    if test -f ${sampleName}_delly.bcf; then
        echo "${sampleName}_delly.bcf exist"
    else
        cp ${baseDir}/bin/delly/mock_output_delly.bcf ${sampleName}_delly.bcf
    fi

    bcftools view ${sampleName}_delly.bcf > ${sampleName}_delly.vcf
    bgzip ${sampleName}_delly.vcf && tabix -p vcf ${sampleName}_delly.vcf.gz
    """
}

/*******************************************************
* Part 4: Filter SNPs and indels from VCF and perform QC
********************************************************/

// Process 4A: filter SNPs from the VCF and apply QC filters (GATK recommendations)
process '4A_VCF_filtering_SNPs' {
tag '4A'
errorStrategy error_strategy
maxRetries max_retries
conda 'bioconda::bcftools=1.8'
publishDir outDir + '/filtered', mode: 'copy'
 input:
    file genome from genome_4A
    file genome_idx from genome_idx_4A
    file dict from fasta_dict_4A
    file vcf from GATK_4A
 output:
    file "${sampleName}_filtered_snps.vcf" into snps_vcf_filtered_4C, snps_vcf_filtered_4D
    file "${sampleName}_flagged_snps.vcf"
    file ".command.*"
 script:
     """
     # extract SNPs only from VCF
     gatk SelectVariants -R ${genome} -V ${vcf} \
     --select-type-to-include SNP --output ${sampleName}_snps_only.vcf

     # recommended: QUAL > 30, QD < 2, FS > 60, MQ < 40, MQRankSUM < -12.5, ReadPosRankSum < -8.0
     # GATK variant filtration step, mark bad calls with flags
     gatk VariantFiltration -R $genome -V ${sampleName}_snps_only.vcf --filter-expression \
     "QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0" \
     --filter-name "FAILED" --output ${sampleName}_flagged_snps.vcf

     # filter SNPS on PASS filter, this file will be used for SNP-it
     bcftools view -f .,PASS ${sampleName}_flagged_snps.vcf > ${sampleName}_filtered_snps.vcf
     """
}

//Process 4B: filter the indels from GATK and Delly
process '4B_VCF_filtering_indels' {
    tag '4B'
    errorStrategy error_strategy
    maxRetries max_retries
    conda 'bioconda::bcftools=1.8'
    publishDir outDir + '/filtered', mode: 'copy'
    input:
        file genome from genome_4B
        file genome_idx from genome_idx_4B
        file dict from fasta_dict_4B
        file gatk_vcf from GATK_4B
        file delly_vcf from delly_vcf
    output:
        file "${sampleName}_gatk_delly_merged.sorted.vcf" into indels_vcf_filtered_4C
        file "${gatk_vcf.baseName}_flagged_indels.vcf"
        file "${delly_vcf.simpleName}_filt.vcf"
        file ".command.*"
    script:
    """
    # extract indels only from vcf
    gatk SelectVariants -R $genome -V ${gatk_vcf} \
    --select-type-to-include INDEL --output ${gatk_vcf.baseName}_indels_only.vcf

    # Variant filtration step for GATK VCF, recommended: QUAL> 30, QD< 2 FS> 60, MQ< 40
    # GATK variant filtration step, mark bad calls with flags
    gatk VariantFiltration -R $genome -V ${gatk_vcf.baseName}_indels_only.vcf --filter-expression \
    "QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "FAILED" --output ${gatk_vcf.baseName}_flagged_indels.vcf

    # second, throw away the failed indels
    bcftools view -f .,PASS ${gatk_vcf.baseName}_flagged_indels.vcf > ${gatk_vcf.baseName}_filt_indels.vcf

    # Variant filtration step for Delly VCF, first only keep 1/1 Genotypes
    bcftools filter -i 'GT="1/1"' ${delly_vcf} > ${delly_vcf.simpleName}_gt_filt.vcf

    # Second, throw away low qual
    bcftools view -f .,PASS ${delly_vcf.simpleName}_gt_filt.vcf > ${delly_vcf.simpleName}_filt.vcf

    # merge DELLY results with GATK results
    ${consensus_vcf} --vcf1 ${gatk_vcf.baseName}_filt_indels.vcf --vcf2 ${delly_vcf.simpleName}_filt.vcf --stdout > ${sampleName}_gatk_delly_merged.vcf

    # sort GATK/DELLY VCF file
    $picard SortVcf I=${sampleName}_gatk_delly_merged.vcf O=${sampleName}_gatk_delly_merged.sorted.vcf
    """
    }

//Process 4C: merge the SNP and indel VCF files from GATK and Delly
process '4C_merge_snps_indels' {
    tag '4C'
    errorStrategy error_strategy
    maxRetries max_retries
    conda 'bioconda::samtools=1.9 bioconda::bcftools=1.8'
    publishDir outDir + '/merged', mode: 'copy'
    input:
        file indels from indels_vcf_filtered_4C
        file snps from snps_vcf_filtered_4C
    output:
        file "${sampleName}_snps_indels_filt.vcf" into filtered_vcf_5A
        file ".command.*"
    script:
    """

    # compress and index indels vcf (required for merging)
    bgzip ${indels} && tabix -p vcf ${indels}.gz

    # compress and index snps vcf (required for merging)
    bgzip ${snps} && tabix -p vcf ${snps}.gz

    # merge two vcf files together, select only the variants with the 'PASS' filter
    bcftools concat -a ${indels}.gz ${snps}.gz -o ${sampleName}_snps_indels_filt.vcf
    """
}

// Process 4D: Create alternative fasta file using reference sequence and PASSED SNPs of strain
process '4D_Create_SNPfasta' {
    tag '4D'
    errorStrategy error_strategy
    maxRetries max_retries
    conda 'bioconda::samtools=1.9'
    publishDir outDir + '/SNPFasta', mode: 'copy'
    input:
        file vcf from snps_vcf_filtered_4D
        file genome from genome_4D
        file genome_idx from genome_idx_4D
        file dict from fasta_dict_4D
    output:
        file "${sampleName}*_altref.fasta.gz" into altfasta_4E
        file ".command.*"
    script:
        """
        # compress and index indels vcf (required for FastaAlternateReferenceMaker)
        bgzip ${vcf} && tabix -p vcf ${vcf}.gz

        # obtain full length genome with SNPs imported
        gatk FastaAlternateReferenceMaker -R ${genome} -V ${vcf}.gz -O ${sampleName}_altref.fasta

        # compress the fasta file
        bgzip ${sampleName}_altref.fasta
        """
}

// Process 4E: Lineage determination for Mycobacteria of the Tuberculosis complex
// will only be performed if accesionnumber == NC_000962.3
process '4E_SNPIT' {
    tag '4E'
    publishDir outDir + '/SNPIT', mode: 'copy'
    input:
        file fasta from altfasta_4E
        file meta from meta_json_snpit
        val acc from acc_4E
    output:
        file "${meta.baseName}.json"
        file ".command.*"
    when:
      "${acc}" == "NC_000962.3"
    script:
        """
        $snpit --input ${fasta} | $snpit_reporter --JSON ${meta} --sample ${sampleName}
        """
}

/**************************************************
 * PART 5: Determine drug resistance of given sample
 **************************************************/

 // Process 5A: annote the VCF for mapping to the resistance database
process '5A_SNPeff_annotation' {
    tag '5A'
    errorStrategy error_strategy
    maxRetries max_retries
    conda 'bioconda::snpeff=4.3.1t'
    publishDir outDir + '/annot', mode: 'copy'
    input:
        file vcf from filtered_vcf_5A
        file gff from gff_5A
        file gbk from gbk_5A
        file genome from genome_5A
        val acc from acc_5A
    output:
        file "${vcf.baseName}_annot.vcf" into GATK_vcf_annot
        file "snpEff_summary.html"
        file ".command.*"
    script:
        """
        snpeff_runtime="${snpeff_base}/data/${acc}/${acc}.config"
        if [ -d ${snpeff_base}/data/${acc} ]
        then
          snpEff -v -interval ${gff} ${acc} ${vcf} -ud 400 -c \${snpeff_runtime} > ${vcf.baseName}_annot.vcf;
        else
          echo "Please build the snpeff database first."
          exit 1
        fi
        """
}

// Process 5B: map the VCF to the resistance database to determine if resistance related mutations are present
process '5B_profile_resistance' {
  tag '5B'
  publishDir outDir + '/mutations', mode: 'copy'
  input:
    file annot from GATK_vcf_annot
    file meta from  meta_json_resistance
    file gff from gff_5B
    val acc from acc_5B
  output:
    file "${meta.baseName}.json" into meta_6A
    file ".command.*"
  script:
    """
    # for extractfields on allele frequency
      $snpsift extractFields ${annot} "ANN[*]" "DP" "GEN[*].AD" | \
      $resis_reporter --acc ${acc} --gff ${gff} --resDB_dir ${resDB} --json $meta
    """
}

// Process 5C: check the depth in the BAM file to determine missed resistance positions
process '5C_check_no_coverage' {
  tag '5C'
  publishDir outDir + '/zerodepth', mode: 'copy'
  conda 'bioconda::samtools=1.20'
  input:
    file bam from bam_5C
    file meta from meta_json_nocov
    val acc from acc_5C
  output:
    file "${sampleName}_zerodepth.tsv"
    file "${meta.baseName}.json"
    file ".command.*"
  when:
      "${acc}" == "NC_000962.3"
  script:
    """
    # Determine depth for all positions on the reference genome
    echo ${bam} > samtools_depth_input.txt
    samtools depth -f samtools_depth_input.txt -aa > ${sampleName}.coverage

    # List relevant resistance positions with zero depth
    $check_depth -s ${sampleName} -c ${sampleName}.coverage -r ${resDB}/${acc}_resistance_db_genomic.tsv -o . -j ${meta}
    """
}

/*********************************************************************************************************
 * PART 6: Generate a report with information about the identified organisms and mutations on target genes
 *********************************************************************************************************/

// Process 6A: store versions of the databases used in the meta json
process '6A_store_versions' {
  tag '6A'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'conda-forge::jq=1.6'
  publishDir outDir + '/database_versions', mode: 'copy'
  input:
    val acc from acc_6A
    file meta from meta_6A
  output:
    file "*.json" into report_json
  script:
    """
    res_version=\$(readlink -v ${resDB}/${acc}_resistance_db.tsv)
    cat ${meta} | \
    jq --arg key "genome_db" --arg value "${genome_db}/${acc}" '.database_versions += {(\$key): \$value}' | \
    jq --arg key "centrifuge_ref" --arg value "${centrifuge_ref}" '.database_versions += {(\$key): \$value}' | \
    jq --arg key "centrifuge_wgs" --arg value "${centrifuge_WGS}" '.database_versions += {(\$key): \$value}' | \
    jq --arg key "centrifuge_cont" --arg value "${centrifuge_cont}" '.database_versions += {(\$key): \$value}' | \
    jq --arg key "resDB" --arg value "\${res_version}" '.database_versions += {(\$key): \$value}' | \
    jq --arg key "hsp65" --arg value "${hsp65}" '.database_versions += {(\$key): \$value}' > \
    ${sampleName}_report.json
    """
}
 
 
// Process 6B: generate a report for interpretation by the clinician (or for research purposes)
process '6B_report' {
  tag '6B'
  errorStrategy error_strategy
  maxRetries max_retries
  conda 'conda-forge::jinja2=2.11.1 conda-forge::weasyprint=51 conda-forge::simplejson=3.17.0 conda-forge::matplotlib=3.1.2 conda-forge::pandas=1.0.1'
  publishDir outDir + '/report', mode: 'copy'
  input:
    file meta from report_json
  output:
    file "${meta.baseName}.json"
    file "${sampleName}.html"
    file "${sampleName}.pdf"
    file ".command.*"
  script:
    """
      $add_stats --JSON $meta
      $reporter --JSON $meta --sample ${sampleName}
    """
}