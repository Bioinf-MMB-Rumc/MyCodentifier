snpeff_base        = "$baseDir/bin/snpEff"
genome_db          = "$baseDir/db/genomes"

// Parsing the input parameters
accession       = "$params.accession"

// databases
template_file    = "$baseDir/bin/snpEff/snpEff.config"
snpeff_base      = "$params.snpeff_base"
conversion_table = "$baseDir/db/centrifuge_db/conversion_table.txt"
genome_dir       = "$params.genome_db/$accession"

// Tools paths and command prefixes
generate_db      = "$baseDir/bin/gen_snpeff_db.py"


log.info """
NEXTFLOW Build databases for MyCodentifier
================================
accession     : $params.accession
snpeff_base    : $params.snpeff_base
================================
"""

/*********************************************************************************
 * PART 1: Building the genome database
 *********************************************************************************/

// Process 1A: Download the fasta, genbank (GFF) files from NCBI if not already present in target folder
process '1A_Build_genome_database' {
    conda 'conda-forge::requests=2.23.0 conda-forge::parallel=20200322'
    tag '1A'
    publishDir genome_dir, mode: 'copy'
    input:
    output:
        file "${accession}_genomic.fna" into genome_2A
        file "${accession}_genomic.gff" into gff_2A, gff_3A
        file "${accession}_genomic.gbff" into gbk_2A
        file "${accession}_filt_genomic.fa"
        file "${accession}_filt_genomic.fa.fai"
        file ".command.*"
    script:
    """

    mkdir -p ${genome_dir}

    # get the genome files from NCBI if not already present in the target directory
    fetch_refs.py --accession ${accession} --table $conversion_table --format all --outdir ${genome_dir}

    # unzip all the genome files in the current directory (in parallel)
    gunzip -c ${genome_dir}/*gff.gz > ${accession}_genomic.gff
    gunzip -c ${genome_dir}/*gbff.gz > ${accession}_genomic.gbff
    gunzip -c ${genome_dir}/*fna.gz > ${accession}_genomic.fna


    #parallel 'zcat {} > {}' ::: ${genome_dir}/*.gz

    # replace type-material INFO (2020-02-18, recently added by NCBI causes problems with snpEff)
    sed -i 's/type-material/typematerial/g' ${accession}_genomic.gff
    sed -i 's/collection-date/collectiondate/g' ${accession}_genomic.gff
    sed -i 's/lat-lon/latlon/g' ${accession}_genomic.gff
    sed -i 's/isolation-source/isolationsource/g' ${accession}_genomic.gff
    sed -i 's/isolation-source/isolationsource/g' ${accession}_genomic.gff
    sed -i 's/culture-collection/culturecollection/g' ${accession}_genomic.gff
    sed -i 's/old-name/oldname/g' ${accession}_genomic.gff
    sed -i 's/collected-by/collectedby/g' ${accession}_genomic.gff
    sed -i 's/nat-host/nathost/g' ${accession}_genomic.gff
    sed -i 's/plasmid-name/plasmidname/g' ${accession}_genomic.gff


    # filter accession code chromosomal DNA (exlude Plasmidic DNA)
    cp ${accession}_genomic.fna ${accession}_filt_genomic.fa
    samtools faidx ${accession}_filt_genomic.fa
    """
}



 /*********************************************************************************
 * PART 2: Building the snpEff database
 *********************************************************************************/

// Process 2A: annote the VCF for mapping to the resistance database
process '2A_build_snpEff_database' {
    tag '2A'
    conda 'bioconda::snpeff=4.3.1t'
    input:
        file gff from gff_2A
        file gbk from gbk_2A
        file genome from genome_2A
    output:
    script:
        """
        cp -r ${snpeff_base} ./

        # run db generation script, if db exist, don't build, else build
        pass=\$(${generate_db} --acc ${accession} --genbank ${gbk} --template ${template_file} \
        --db ${snpeff_base}/data/${accession} --outDir ${snpeff_base}/data/${accession} )
        if [ "\$pass" == "0" ];
        then
        # check if directory exist and makedir, copy genome to destination ,build snpeff db
            [ ! -d ${snpeff_base}/data/${accession} ] && mkdir -p ${snpeff_base}/data/${accession}
            cat ${genome} > ${snpeff_base}/data/${accession}/sequences.fa
            mv ${snpeff_base}/data/${accession}/runtime.config ${snpeff_base}/data/${accession}/${accession}.config
            gzip < ${gff} > ${snpeff_base}/data/${accession}/genes.gff.gz && snpEff build -gff3 -v ${accession} -c ${snpeff_base}/data/${accession}/${accession}.config;
        elif [ "\$pass" != "1" ];
        then
            echo \${pass};
            exit 1;
        fi
        """
}

