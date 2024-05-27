#!/usr/bin/env nextflow

nextflow.enable.dsl = 1

params.help = false

def parameters_expected = ['path', 'ref_path', 'ref_list', 'singleEnd',
                           'outdir',
			   'skip_profile', 'run_minpath', 'max_strain_rank', 'strain_rank',
			   'single-end', // why this is in?
			   'awsregion','awsqueue',
			   'max_memory', 'max_cpus', 'max_time',
			   'pipelineVersion', 'pipeline-version', 'tracedir', 'help'
			   ] as Set


def helpMessage() {
  // adapted from nf-core
    log.info"""

    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run main.nf --path FODLER_FOR_READS --ref_path FOLDER_FOR_REFERENCEs
    Mandatory arguments:
      --path                        Path to a folder containing all input fastq files (this will be recursively searched for *fastq.gz files)
      --ref_path                    Path to PanPhlAn reference databases
    Optional arguments:
      --ref_list                    A list of reference species to profile. Each line is a full name of a species in the database, e.g. Faecalibacterium prausnitzii
      --singleEnd                   The fastq files are single-end reads [Default: false]
    StrainPanDAR options:
      --skip_profile                Skip decomposition algorithm in strainpandar [Default: false]
    
    MinPath options:
      --run_minpath                 Run MinPath [default: false]
    Other:
      --outdir                      The output directory where the results will be saved (Default: panphlan_out)
    AWSBatch:
      --awsregion                   The AWS Region for your AWS Batch job to run on
      --awsqueue                    The AWS queue for your AWS Batch job to run on
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
if (!params.path){
    helpMessage()
    log.info"""
    [Error] --path is required
    """.stripIndent()
    exit 0
}
if (!params.ref_path){
    helpMessage()
    log.info"""
    [Error] --ref_path is required
    """.stripIndent()
    exit 0
}
if (!params.ref_list){
    // use all genomes
    ch_ref = Channel
        .fromPath( params.ref_path + '/**_pangenome.csv' )
        .map{ x -> params.ref_path + '/' + x.name.split("_")[1] }
	.map{ x -> tuple(x.split("/").last(), file(x))}
}else{
    ch_ref = Channel
        .fromPath( params.ref_list )
        .splitText()

        .map{ x -> tuple(params.ref_path + '/' +
	      	   	 x.replaceAll(' ','-').trim(), (params.ref_path =~ /[0-9]+\/*$/)[0].replaceAll("/",""))  }
	.map{ x -> tuple(x[0].split("/").last() + "-" + x[1], file(x[0] + "-" + x[1])) }
}


def parameter_diff = params.keySet() - parameters_expected
if (parameter_diff.size() != 0){
   exit 1, "Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
   }

ch_ref
    .into { ch_ref_map; ch_ref_profile }

Channel
    .fromFilePairs( params.path + '/**{1,2}.f*q*', size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
    .combine(ch_ref_map)
    .set{ ch_panphlan_input }

process panphlanMap {
    tag "${prefix}|${species}"
    //cache 'deep'

    publishDir params.outdir, mode: 'copy'

    input:
    set prefix, file(reads), species, file(ref_path) from ch_panphlan_input

    output:
    set species, file("${species}_mapping/${prefix}_${species}*") into ch_panphlanmap

    """
    ext=`echo $reads | sed 's/.*\\.//'`
    concat=\$([ "\$ext"  == "gz" ] && echo "zcat" || echo "cat")
    \$concat $reads \\
    | panphlan_map.py \\
        --i_bowtie2_indexes ${ref_path} \\
        -c ${species} \\
        -o ${species}_mapping/${prefix}_${species}.csv \\
        -p $task.cpus
    """
}


ch_panphlanmap
    .groupTuple()
    .set {ch_mapped}

process mergeProfile {
    tag "${prefix}"

    publishDir params.outdir, mode: 'copy'

    input:
    set prefix, file(mapped) from ch_mapped

    output:
    set prefix, file("${prefix}.counts.csv") into ch_merged

    """
    Rscript /opt/bin/merge_tables.r -p "*bz2" -o ${prefix}.counts.csv
    """
}

ch_merged
    .join(ch_ref_profile)
    .set { ch_species_counts }

process runStrainPanDAR {
    tag "${prefix}"

    publishDir "${params.outdir}/${prefix}_strainpandar_out", mode: 'copy', pattern: '*{pdf,rds,csv}'

    when:
    params.skip_profile == false

    input:
    set prefix, file(counts), file(pangenome) from ch_species_counts

    output:
    set prefix, file("${prefix}.strainpanda.*pdf"), file("${prefix}.strainpanda.*rds"), file("${prefix}.strainpanda.*csv") into ch_strainpandar
    file("${prefix}.*ko.txt") into ch_ko
    file("*.{pdf,rds,csv,txt}") into ch_all

    """
    Rscript /opt/bin/run_strainpandar.r -c $counts -r $pangenome -o ${prefix}.strainpanda -t $task.cpus -m $params.max_strain_rank -n $params.strain_rank
    """
}


process runMinPath {
    tag "${prefix}"

    publishDir "${params.outdir}/${prefix}_minpath_out", mode: 'copy'

    when:
    params.skip_profile == false && params.run_minpath == true

    input:
    set prefix, file(ko) from ch_ko.flatten().map {file -> tuple( file.name.split("\\.")[0], file)}

    output:
    set prefix, file("*minpath") into ch_minpath

    """
    python /MinPath/MinPath1.4.py -ko $ko -report ${ko}.tmp
    grep "minpath 1" ${ko}.tmp | cut -f21- -d" " | awk -v name=$ko '{print name"\t"\$0}' > ${ko}.minpath
    """

}


