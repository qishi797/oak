#!/usr/bin/env nextflow


// params. use --parameter_name to change parameter
params.paired_reads = './data/reads/R10/*{1,2}.fq.gz' // remember to change this to null. Use example --paired_reads='./data/reads/*{1,2}.fq.gz'
params.reference_genome= "${baseDir}/data/references/Hazelnut_CavTom2PMs-1.0/GCF_932294415.1" //This path should be the full path to reference genome. 
params.reference_name = "GCF_932294415.1_dhQueRobu3.1_genomic.fa"
params.pipeline_loc = "/rds/projects/h/hejq-msc-2024-wgbs/Test_Sbatch/" //full path such as full/path/to/Hazexplore/
params.results = "./results"
params.temps = "${baseDir}/temps"
params.index_requirement = 0 //change this to null 
params.parallelize = 4
params.threads = 6

// Help function with useful parameters.
def help = params.help ?: false
if (params.help){
    log.info """\
    
    Command options:
    
    --paired_reads=<pattern>      Path to paired-end reads in FASTQ format.  Use the pattern "./data/reads/*{1,2}.fq.gz". By default reads are looked into the ./data/reads/ directory
    --reference_genome=<path>     Full path to the reference genome in FASTA format.
    --reference_name=<name>       Name for the reference genome. Make sure it contains the .fa reference file. 
    --results=<directory>         Directory to store the pipeline results (default: ./results).
    --index_requirement=<value>   Specify an integer (0 or 1) to indicate if indexing the reference genome is required (0: not required, 1: required).
    --parallelize=<value>         Specify the level of parallelization (default: 1).
    --threads=<value>             Specify the number of threads to use for parallel tasks (default: 6).
    --help                        Display this help message and exit.
    """
    exit (" ")
}


// Global variables
trimmed_outputs= "${baseDir}/data/trimmed/"
indexed_reference = "${baseDir}/data/references/"
temps = "${baseDir}/temps/"

// Checks if fundamental parameters have been specified.
if (params.paired_reads == null) {
      println("Specify reads using paired_reads")
      exit(1)
}

if (params.index_requirement == null){
    println("Specify integer for indexing the reference genome is necessary. (0: non required, 1: required)")
    exit(1)
}

// Info for the user.
log.info """\

    HAZEXPLORER - NF PIPELINE
    ===================================
    
    Reads: ${params.paired_reads}
    Reference genome: ${params.reference_genome}

    """


// Trims the reads
process TRIM {
    tag {sampleId}
    publishDir "${trimmed_outputs}/${sampleId}"
    
    input:
    tuple val(sampleId) , path(reads) 

    output:
    tuple val(sampleId), path ("*.fastq.gz")

    script:
    def (read1, read2) = reads
    """
    set -e

    module purge; module load bluebear
    module load bear-apps/2021b/live
    module load fastp/0.23.2-GCC-11.2.0

    mkdir -p "${trimmed_outputs}${sampleId}"

    # Trim paired sequences using fastp (approximately 45 seconds)
    fastp -i "${read1}" -I "${read2}" \
    -o "${sampleId}_trimmed_1.fastq.gz" \
    -O "${sampleId}_trimmed_2.fastq.gz"
    """
}


//Indexes reference genome based on bowtie2 algorithm.
process INDEX{
    
    publishDir "${indexed_reference}/${reference_name}/"
    
    input: 
        path reference_genome
    output: 
        path "${params.reference_name}/*"
    script:

    """
    set -e

    module purge 
    module load bluebear 
    module load bear-apps/2021b
    module load Bismark/0.24.2-foss-2021b

    
    bismark_genome_preparation --bowtie2 --verbose ${reference_genome} / 
    """
// change path to aligner so it not /usr/bin/bowtie2, but just bowtie2.  
}

// creates sequence dictionary for reference genome for use during BIS SNP
process PICARD_DICT{

    publishDir "${params.reference_genome}/"

    input: 
        path "${params.reference_name}"
        
    //output
        //path "${params.reference_genome}/"

    script:


    """
    set -e 

    module purge; module load bluebear
    module load bear-apps/2022a/live
    module load bear-apps/2022b/live
    module load Java/17.0.6
    
    if [! -f ${params.reference_name}.dict]; then
    
        java -Xmx4g -jar ${params.pipeline_loc}/tools/picard.jar CreateSequenceDictionary \
        R=${params.reference_genome}/${params.reference_name} \
        TRUNCATE_NAMES_AT_WHITESPACE=true NUM_SEQUENCES=2147483647 VERBOSITY=INFO QUIET=false \
        VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false \
        CREATE_MD5_FILE=false 

    fi
    """
}

process FAST_QC{
    tag {sampleId}
    publishDir "${params.results}/QC/${sampleId}"
    
    input:
    tuple val(sampleId) , path(reads) 

    output: 
    path "*.html"
    path "*.html"


    script:
    def (read1, read2) = reads
    """
    set -e

    module purge; module load bluebear
    module load bear-apps/2021b/live
    module load FastQC/0.11.9-Java-11
    
    
    mkdir -p "${params.results}/QC/${sampleId}"

    fastqc ${read1} --threads ${params.threads} --quiet true --output ${sampleId}_QC_1
    fastqc ${read2} --threads ${params.threads} --quiet true --output ${sampleId}_QC_2
    """
}

//Aligns reads using bismark and bowtie2 
process ALIGNMENT {
    tag {sampleId}
    
    input:
    tuple val(sampleId) , path(reads) 
    path indexed_reference_directory 
    
    output:
    tuple val(sampleId), path ("*.bam"), path ("*.txt")
    

    script:
    def (trimmedRead1, trimmedRead2) = reads
    """
    set -e

    module purge; module load bluebear 
    module load bear-apps/2021b
    module load Bismark/0.24.2-foss-2021b

    mkdir -p "${params.temps}/Alignments/${sampleId}"

    bismark --bowtie2 -p ${params.threads} --multicore ${params.parallelize} --genome ${indexed_reference_directory} \
     -1 ${trimmedRead1} -2 ${trimmedRead2} --bam 
    """
}


//Create Picard Sorting in preparation for SNP calling

process PICARD{
    tag {sampleId}


    input:
    tuple val(sampleId), path (bam_file_in), path (bismark_report)

    output: 
    tuple val(sampleId), path ("*.bam")
    
    script: 
    def bam_file = bam_file_in
    """
    set -e 

    module purge; module load bluebear
    module load bear-apps/2022a/live
    module load bear-apps/2022b/live
    module load Java/17.0.6


    java -Xmx4g -jar  ${params.pipeline_loc}/tools/picard.jar AddOrReplaceReadGroups \
    I=${bam_file} \
    O=${sampleId}_pic_uns.bam \
    RGID=${sampleId}_RG \
    RGLB=Unknown \
    RGPL=ILLUMINA \
    RGPU=Unknown \
    RGSM=${sampleId} \
    CREATE_INDEX=false

    """

} 



//carries out sorting of BAM files edited with PICARD. 
process SORTING{
    tag {sampleId}

    publishDir "${params.results}/Alignments/${sampleId}/"

    input:
    tuple val(sampleId), path (unsorted_bam_file)

    output: 
    tuple val(sampleId), path ("*.bam"), path ("*.bai")
    
    script: 
    def bam_file = unsorted_bam_file
    """
    set -e

    module purge; module load bluebear
    module load bear-apps/2022a/live
    module load bear-apps/2022b/live
    module load Java/17.0.6


    java -Xmx4g -jar  ${params.pipeline_loc}/tools/picard.jar SortSam \
    I=${bam_file} \
    O=${sampleId}_sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=${params.pipeline_loc}/temps/

    """


}


//Create a process for SNP calling using BisSNP.
process BIS_SNP {
    tag { sampleId }
    publishDir "${params.results}/results/${sampleId}/"

    input:
    tuple val(sampleId), path(sorted_bam_file), path(bai_file)

    output:
    tuple val(sampleId), file("${sampleId}*.vcf"), file("${sampleId}*.txt")

    script:
    def bam_file = sorted_bam_file
    def bai_index = bai_file
    """
    set -e

    module purge; module load bluebear
    module load bear-apps/2022b/live
    module load Java/17.0.6
    module load SAMtools/1.17-GCC-12.2.0
    module load GATK/4.4.0.0-GCCcore-12.2.0-Java-17

    # Calls SNPs using BisSNP
    java -Xmx4g -jar ${params.pipeline_loc}/tools/BisSNP-0.90.jar -R ${params.reference_genome}/*.fa \
    -nt 10 -T BisulfiteGenotyper -I ${bam_file} \
    -vfn1 ${sampleId}_cpg.raw.vcf -vfn2 ${sampleId}_snp.raw.vcf

    # Add command that Generates a summary table and graph for SNP amount found at each chromosome.
    
    """
}



//Carries out prep for  CGMapTools used in bayesian mode.
process  CGMAP_PREP {
    tag { sampleId }
    publishDir "${params.results}/results/${sampleId}/"

    input:
    tuple val(sampleId), path(sorted_bam_file), path(bai_file)
    path reference_genome

    output:
    tuple val(sampleId), path ("*.ATCGmap.gz"), path ("*.CGmap.gz")
    
    script: 
    def bam_file = sorted_bam_file
    def bai_index = bai_file
    def reference_genome = reference_genome

    """
    set -e

    module purge; module load bluebear
    module load bear-apps/2022a/
    module load CGmapTools/0.1.3-foss-2022a


    #Convert BAM file into input files for CGmap tools
    cgmaptools convert bam2cgmap -b ${bam_file} -g ${reference_genome} -o ${sampleId}

    """


}


//Create a process for SNP calling using CGMapTools in bayesian mode.
process CGMAP_TOOLS {
    tag { sampleId }
    publishDir "${params.results}/results/${sampleId}/"

    input:
    tuple val(sampleId), path(ATCGmap_file), path(CGmap_file)
    

    output:
    tuple val(sampleId), path("*.vcf"), path("*.snv")

    script:
    def ATCGmap_file = ATCGmap_file
    """
    set -e

    module purge; module load bluebear
    module load bear-apps/2022a/
    module load CGmapTools/0.1.3-foss-2022a

    # Carrys out SNP analysis using bayesian method. Slower but more precise.
    cgmaptools snv -i "${ATCGmap_file}" -m bayes -v "${sampleId}_bayes.vcf" -o "${sampleId}_bayes.snv" --bayes-dynamicP

    """

}


// Acts as the MAIN function, running each process in the most optimal way.
workflow{

    //reads mates from a specified directory.
    paired_reads_ch= Channel.fromFilePairs(params.paired_reads, checkIfExists: true)
    paired_reads_ch.view()
    
    paired_trimmed = TRIM(paired_reads_ch)
    FAST_QC(paired_trimmed) //produces a QC report of trimmed reads. 
    PICARD_DICT(params.reference_genome)
   
    //called if reference genome is default or does not need indexing.
    if (params.index_requirement == 0){
        aligned_bam = ALIGNMENT(paired_trimmed, params.reference_genome)
        picard_out = PICARD(aligned_bam)
        sort_out = SORTING(picard_out)
        bis_snp_out = BIS_SNP(sort_out)
        
        CGmaptool_prp = CGMAP_PREP(sort_out, "${params.reference_genome}/${params.reference_name}")
        CGmap_out = CGMAP_TOOLS(CGmaptool_prp)
            
    }
    //called if reference genome is custom and needs to be indexed.
    if (params.index_requirement == 1){
        indexed_reference = INDEX(params.reference_genome)
        aligned_bam = ALIGNMENT(paired_trimmed, params.reference_genome)
        picard_out = PICARD(aligned_bam)
        sort_out = SORTING(picard_out)
        bis_snp_out = BIS_SNP(sort_out)

        CGmaptool_prp = CGMAP_PREP(sort_out, "${params.reference_genome}/${params.reference_name}")
        CGmap_out = CGMAP_TOOLS(CGmaptool_prp)

    }

}