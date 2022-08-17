#!/usr/bin/env nextflow

/*
 * =======================================================================
 * 
 *		Bidirectional Transcription Detection Pipeline
 * 		
 * =======================================================================
 * 
 * This Source Code is distributed under the GPL-3.0 License
 */


/* 
 * 'Bidirectional-Flow' - A Nextflow pipeline for detecting regions of transcription
 * from Nascent RNA sequencing  
 *
 * The pipeline identifies regions of transcription with FStitch, TFit and/or dREG
 * 
 * =============
 * Authors
 * =============
 * Rutendo F. Sigauke : rutendo.sigauke@cuanschutz.edu
 * Lynn Sanford : lynn.sanford@colorado.edu
 */

def helpMessage() {
    log.info"""
    =========================================
     BidirectionalFlow v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile slurm --crams '/project/*.sorted.cram' --workdir '/project/tempfiles' --outdir '/project/'
    Required arguments:
         -profile                      Configuration profile to use. <genome_user>
         --crams                       Directory pattern for cram files: /project/*.sorted.cram (Required if --bams not specified).
         --bams                        Directory pattern for bam files: /project/*.sorted.bam (Required if --crams not specified).
         --workdir                     Nextflow working directory where all intermediate files are saved.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --r1_five_prime		       If input file is paired, specifies if read 1 has the 5 prime end (default R2 is five prime, must be manually determined)

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savestats                    Saves tfit/dreg bidirectional statistics for all samples (default true)
        --savebam                      Saves sorted bam files (if cram as input)
        --savetfitbam                  Saves multimapped read filtered bamfiles used for tfit
        --savebg                       Saves bedgraph files for tfit and/or dreg (whichever is run)
        --savebw                       Saves dreg bigwig files
        --savebidirs                   Saves bedfiles for promoter/intronic/intergenic bidirectional subsets

    Analysis Options:
        --gene_count                   Run featureCounts to obtain stranded and unstranded gene counts over an annotation.
        --bidir_count                  Run featureCounts to obtain stranded and unstranded bidirectional counts over an annotation (SAF file)
        --fstitch                      Run FStitch. If used, you must also specify FS_path and FS_train params.
        --tfit                         Run Tfit bidir. If used, you must also specify the Tfit_path parameter.
        --tfit_prelim		       Run Tfit bidir. If used, you must also specify the Tfit_path parameter. Not compatible with --prelim_files flag.
        --tfit_split_model             Run Tfit model separately on <5kb prelim regions and 5-10kb prelim regions (with different k parameters). Not compatible with --tfit_model or --tfit flags, or --prelim_process=false
        --tfit_model		       Run Tfit full model. If used, must specify the Tfit path parameter AND have prelim files from --tfit_prelim process or previous run via the --prelim_files flag. Not compatible with --tfit flag.
        --prelim_files		       Directory pattern for tfit prelim files: /project/*-1_prelim_bidir_hits.bed (required for --tfit_model if --tfit_prelim is not also specified)
        --prelim_process               Run the prelim processing step for Tfit (default true)
        --dreg                         Produce bigwigs formatted for input to dREG.
        --dreg_results                 Do coverage filtering on existing dreg results
    """.stripIndent()
}

// Configure Variables

params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.tfit_prelim_run = "$baseDir/bin/tfit_prelim.sh"
params.tfit_model_run = "$baseDir/bin/tfit_model.sh"
params.prelim_filter = "$baseDir/bin/prelim_filter.py"
software_versions = Channel.create()

import java.text.SimpleDateFormat
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd")
output_date =  sdf.format(date)

String output_date = new java.text.SimpleDateFormat("yyMMdd").format(new Date())

// Header log info
log.info """=======================================================
Bidirectional-Flow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'Bidirectional-Flow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = workflow.runName
if(params.crams) summary['Crams']            = params.crams
if(params.bams) summary['Bams']              = params.bams
summary['Genome Ref']       = params.genome
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Five prime read']  = params.r1_five_prime ? 'Read 1' : 'Read 2 or single end'
summary['Output dir']       = params.outdir
summary['FStitch']          = params.fstitch ? 'YES' : 'NO'
summary['Tfit']             = params.tfit ? 'YES' : 'NO'
summary['Tfit prelim']      = params.tfit_prelim ? 'YES' : 'NO'
summary['Tfit model']       = params.tfit_model ? 'YES' : 'NO'
summary['Tfit split model'] = params.tfit_split_model ? 'YES' : 'NO'
summary['Tfit prelim files']	= params.prelim_files ? 'YES' : 'NO'
summary['Process Tfit prelims'] = params.prelim_process ? 'YES' : 'NO'
summary['dREG']             = params.dreg ? 'YES' : 'NO'
summary['dREG postprocess'] = params.dreg ? 'YES' : params.dreg_results ? 'YES' : 'NO'
summary['Gene counting']    = params.gene_count ? 'YES' : 'NO'
summary['Bidir counting']   = params.bidir_count ? 'YES' : 'NO'
if(params.fstitch)summary['FStitch dir']      = params.fstitch_path
if(params.fstitch)summary['FStitch train']    = params.fstitch_train
if(params.tfit)summary['Tfit dir']      = params.tfit_path
summary['Save bidir stats'] = params.savestats ? 'YES' : 'NO'
summary['Save bamfiles']    = params.savebam ? 'YES' : 'NO'
summary['Save Tfit bamfiles']    = params.savetfitbam ? 'YES' : 'NO'
summary['Save bedgraphs']   = params.savebg ? 'YES' : 'NO'
summary['Save bigwigs']     = params.savebw ? 'YES' : 'NO'
summary['Save bidir subsets']   = params.savebidirs ? 'YES' : 'NO'
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="


println "\nSTARTING PIPELINE"

// PART 0: Outputting software versions

println "[Log 0]: Getting software versions"

process get_software_versions {
    time '1h'

    output:
    stdout into software_versions

    script:
    """
    printf "bidirectionalflow_version: %s\n" ${params.version}
    printf "nextflow_version: %s\n" ${workflow.nextflow.version}
    printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}')
    printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}')
    printf "openmpi_version: %s\n" \$(ompi_info | head -2 | tail -1 | awk '{print \$NF}')
    printf "gcc_version: %s\n" \$(gcc --version | head -1 | awk '{print \$NF}')
    printf "fstitch_version: %s\n" \$(${params.fstitch_path} train --version | head -1)
    printf "tfit_version: %s\n" \$(${params.tfit_path} model --version | head -1)
    printf "r_version: %s\n" \$(R --version | head -1 | awk '{print \$3}')
    printf "rsubread_version: %s\n" \$(Rscript -e 'library("Rsubread");packageVersion("Rsubread")' 2>&1 | tail -1 | awk '{print \$NF}')
    printf "boost_version: %s\n" \$(ls -d /Users/\$USER/.local/boost* | head -1 | awk -F "_" '{print \$(NF-2)"."\$(NF-1)"."\$(NF)}')
    printf "dreg_version: %s\n" \$(Rscript -e 'library("dREG");packageVersion("dREG")' 2>&1 | tail -1 | awk '{print \$NF}')
    printf "pipeline_hash: %s\n" ${workflow.scriptId}
    """
}

software_versions.collectFile(name: "software_versions_bidir_${output_date}_${workflow.runName}.yaml", storeDir: "${params.outdir}/pipeline_info")

println "[Log 0]: Software versions complete"


// PART 1: Converting cram files to bam files

if (params.crams) {
  println "[Log 1]: Converting CRAM files to BAM files"
  println "[Log 1]: Genome file being used ..... $params.genome "
  println "[Log 1]: Cram file directory ........ $params.crams"
  println "[Log 1]: Working directory ... $params.workdir"
  println "[Log 1]: Output directory ... $params.outdir"
  cramfiles = Channel
                  .fromPath(params.crams)
                  .map { file -> tuple((file.simpleName + '.sorted'), file)}

  process cram_to_bam {
     cpus 16
     queue 'short'
     memory '5 GB'
     time '1h30m'
     tag "$prefix"

     publishDir "${params.outdir}" , mode: 'copy',
     saveAs: {filename ->
              if (params.savebam && (filename == "${prefix}.sorted.bam"))    "bams/${prefix}.sorted.bam"
              else null
             }

     input:
     tuple val(prefix),file(cram) from cramfiles

     output:
     tuple val(prefix), file("${prefix}.sorted.bam"), file("${prefix}.sorted.bam.bai") into sorted_bam_file, bam_for_dreg, bam_for_gene_counting, bam_for_bidir_counting

     script:
     """
     samtools view -@ 16 -b -1 -T ${params.genome} ${cram} > ${prefix}.sorted.bam
     samtools index ${prefix}.sorted.bam ${prefix}.sorted.bam.bai
     """
  }
} else {

  sorted_bam_file = Channel
                  .fromPath(params.bams)
                  .map { file -> tuple((file.simpleName + '.sorted'), file, (file + '.bai'))}

  bam_for_dreg = Channel
                  .fromPath(params.bams)
                  .map { file -> tuple((file.simpleName + '.sorted'), file, (file + '.bai'))}

  bam_for_gene_counting = Channel
                  .fromPath(params.bams)
                  .map { file -> tuple((file.simpleName + '.sorted'), file, (file + '.bai'))}

  bam_for_bidir_counting = Channel
                  .fromPath(params.bams)
                  .map { file -> tuple((file.simpleName + '.sorted'), file, (file + '.bai'))}

}

process bam_conversion_tfit {
   cpus 16
   queue 'short'
   memory '5 GB'
   time '2h'
   tag "$prefix"

   publishDir "${params.outdir}" , mode: 'copy',
   saveAs: {filename ->
            if (params.savetfitbam && (filename == "${prefix}.mmfilt.sorted.bam"))    "bams/${prefix}.mmfilt.sorted.bam"
            else null
           }

   when:
   params.tfit || params.fstitch || params.tfit_model || params.tfit_prelim || params.tfit_split_model

   input:
   tuple val(prefix), file(bam), file(index) from sorted_bam_file

   output:
   tuple val(prefix), file("${prefix}.mmfilt.sorted.bam"), file("${prefix}.mmfilt.sorted.bam.bai") into bam_for_tfit

   script:
   """
   samtools view -@ 16 -h -q 1 ${bam} | \
       grep -P '(NH:i:1|^@)' | \
       samtools view -h -b > ${prefix}.mmfilt.sorted.bam
   samtools index ${prefix}.mmfilt.sorted.bam ${prefix}.mmfilt.sorted.bam.bai
   """

//   script:
//   """
//   samtools view -@ 16 -h -q 1 ${bam} | \
//       grep -v 'XS:i' | \
//       samtools view -h -b > ${prefix}.mmfilt.sorted.bam
//   samtools index ${prefix}.mmfilt.sorted.bam ${prefix}.mmfilt.sorted.bam.bai
//   """

}

println "[Log 1]: Bam files are ready\n"


// PART 2: Generate bedgraphs

process bedgraphs {
    println "[Log 2]: Generating BEDGRAPHS for TFit and FStitch"
    println "[Log 2]: Genome information ..... $params.genome "
    println "[Log 2]: Chromosome Sizes ....... $params.chrom_sizes"

    tag "$prefix"
    memory '40 GB'
    queue 'short'
    time '4h'
    
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (params.savebg && (filename == "${prefix}.bedGraph"))    "bedgraphs/${prefix}.bedGraph"
             else null
            }

    when:
    params.tfit || params.fstitch || params.tfit_model || params.tfit_prelim || params.tfit_split_model

    input:
    tuple val(prefix), file(bam_file), file(index) from bam_for_tfit

    output:
    tuple val(prefix), file("${prefix}.bedGraph"), file("${prefix}.pos.bedGraph"), file("${prefix}.neg.bedGraph") into fstitch_bg
    tuple val(prefix), file("${prefix}.bedGraph") into prelimtfit_bg, prelimtfit_process_bg, modeltfit_bg, modeltfit_bg_split_max5kb, modeltfit_bg_split_max10kb, post_tfit_bg_split, nqc_bg

    script:
    if (params.singleEnd) {
    """
    genomeCoverageBed \
        -bg \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${bam_file} \
        > ${prefix}.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${bam_file} \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.neg.bedGraph
    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

    """
    } else {
    """
    samtools view \
        -h -b -f 0x0040 \
        ${bam_file} \
        > ${prefix}.first_pair.bam

    samtools view \
        -h -b -f 0x0080 \
        ${bam_file} \
        > ${prefix}.second_pair.bam

    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        > ${prefix}.first_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.first_pair.bam \
        | sortBed \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        > ${prefix}.first_pair.neg.bedGraph

    genomeCoverageBed \
        -bg \
        -split \
        -strand + \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | sortBed \
        > ${prefix}.second_pair.pos.bedGraph
    genomeCoverageBed \
        -bg \
        -split \
        -strand - \
        -g ${params.chrom_sizes} \
        -ibam ${prefix}.second_pair.bam \
        | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' \
        | sortBed \
        > ${prefix}.second_pair.neg.bedGraph

    unionBedGraphs \
        -i ${prefix}.first_pair.pos.bedGraph ${prefix}.second_pair.pos.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.pos.bedGraph

    unionBedGraphs \
        -i ${prefix}.first_pair.neg.bedGraph ${prefix}.second_pair.neg.bedGraph \
        | awk -F '\t' {'print \$1"\t"\$2"\t"\$3"\t"(\$4+\$5)'} \
        > ${prefix}.neg.bedGraph

    cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

    sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

    """
    }
 }


println "[Log 2]: Bedgraph files have been generated\n"


// PART 3: Running FStitch

process FStitch {
    println "[Log 3]: Running FStitch"
    println "[Log 3]: FStitch training file .. $params.fstitch_train"
    println "[Log 3]: FStich source code ..... $params.fstitch_path"

    tag "$prefix"
    memory '50 GB'
    queue 'short'
    time '4h'

    publishDir "${params.outdir}/fstitch/", mode: 'copy', pattern: "*.hmminfo"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "${prefix}.fstitch_seg.ON.merged.bed"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "${prefix}.fstitch_seg.bed"
    publishDir "${params.outdir}/fstitch/segment/", mode: 'copy', pattern: "*.fstitch_seg.{pos,neg}.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "${prefix}.fstitch_bidir.bed"
    publishDir "${params.outdir}/fstitch/bidirs/", mode: 'copy', pattern: "*fstitch_bidir.{short,long}.bed"
    publishDir "${params.outdir}/fstitch/bidirs/hist/", mode: 'copy', pattern: "*.html"
    publishDir "${params.outdir}/fstitch/bidirs/stats/", mode: 'copy', pattern: "*.txt"

    when:
    params.fstitch

    input:
    tuple val(prefix), file(bg), file(pos_bg), file(neg_bg) from fstitch_bg

    output:
    file ("*.hmminfo") into fs_train_out
    tuple val(prefix), file ("*.fstitch_seg.ON.merged.bed") into fs_seg_out
    tuple val(prefix), file ("*.fstitch_bidir.bed") into fs_bidir_out
    file ("*fstitch_bidir.{short,long}.bed") into fs_bidir_short_long_out
    file ("*.html") into fs_bidir_plot_out
    file ("*.txt") into fs_bidir_stats_out

    script:
    """
    ${params.fstitch_path} train \
        -s + \
        -b ${bg} \
        -t ${params.fstitch_train} \
        -o ${prefix}.fstitch.hmminfo

    ${params.fstitch_path} segment \
        -s + \
        -b ${pos_bg} \
        -p ${prefix}.fstitch.hmminfo \
        -o ${prefix}.fstitch_seg.pos.bed

    ${params.fstitch_path} segment \
        -s - \
        -b ${neg_bg} \
        -p ${prefix}.fstitch.hmminfo \
        -o ${prefix}.fstitch_seg.neg.bed

    cat ${prefix}.fstitch_seg.pos.bed \
        ${prefix}.fstitch_seg.neg.bed \
        | sortBed > ${prefix}.fstitch_seg.bed
	    
    cat ${prefix}.fstitch_seg.bed | \
    	grep ON | \
	bedtools merge -i stdin \
	> ${prefix}.fstitch_seg.ON.merged.bed
    
    bidir \
        -b ${prefix}.fstitch_seg.bed \
        -g ${params.genome_refseq} \
        -o ${prefix}.fstitch_bidir.bed \
        -p \
        -s
    """
}

println "[Log 3]: Done Running FStitch\n"


// PART 4: Running Tfit
//Load in external prelim files

if (params.prelim_files) {
    tfit_prelim_out = Channel
        .fromPath(params.prelim_files)
        .map { file -> tuple((file.simpleName + '.sorted'), file)}

// Or run Tfit bidir to make prelim files
} else {

    process tfit_prelim {
        println "[Log 4a]: Running Tfit prelim"

        tag "$prefix"
        memory '100 GB'
        time '6h'
        queue 'short'
        clusterOptions = '-N 1 -c 1'

        publishDir "${params.outdir}/tfit/prelim_logs", mode: 'copy', pattern: "*{log}"
        publishDir "${params.outdir}/tfit/prelim", mode: 'copy', pattern: "*_prelim_bidir_hits.bed"

        when:
        params.tfit_prelim || params.tfit || params.tfit_model || params.tfit_split_model

        input:
        tuple val(prefix), file(bg) from prelimtfit_bg

        output:
        tuple val(prefix), file ("*.sorted-1_prelim_bidir_hits.bed") into tfit_prelim_out
        file ("*.log") into prelimtfit_logs_out

        script:
        """
        ${params.tfit_prelim_run} -t ${params.tfit_path} \
                           -c ${params.tfit_config} \
                           -b ${bg} \
                           -p ${prefix}
        """

    }

println "[Log 4a]: Done Running Tfit prelim\n"

}

//Combine files with same prefix
if (params.tfit_split_model) {
    tfit_prelim_bg_out_preprocess = prelimtfit_process_bg
        .join(tfit_prelim_out)

    process tfit_prelim_process_split {
        println "[Log 4b]: Processing prelim file for split model run"

        tag "$prefix"
        memory '20 GB'
        time '1h'
        queue 'short'

        publishDir "${params.outdir}/tfit/prelim", mode: 'copy', pattern: "*_prelim_coverage_filtered_diced.bed"

        when:
        params.tfit_split_model

        input:
        tuple val(prefix), file(bg), file(prelim) from tfit_prelim_bg_out_preprocess

        output:
        tuple val(prefix), file(bg), file ("*sorted_prelim_coverage_filtered_diced.bed") into tfit_prelim_bg_out
        tuple val(prefix), file(bg), file ("*5kb_prelim_coverage_filtered_diced.bed") into tfit_prelim_bg_out_split_max5kb
        tuple val(prefix), file(bg), file ("*10kb_prelim_coverage_filtered_diced.bed") into tfit_prelim_bg_out_split_max10kb

        script:
        """
        python3 ${params.prelim_filter} \
                -p ${prelim} \
                -b ${bg} \
                -s ${prefix} \
                -o . \
                -g ${params.filtered_refseq} \
                -c ${params.chrom_sizes}

        awk '{if (\$3-\$2 < 5000) print \$0}' \
            ${prefix}_prelim_coverage_filtered_diced.bed \
            > ${prefix}_max5kb_prelim_coverage_filtered_diced.bed

        awk '{if (\$3-\$2 >= 5000) print \$0}' \
            ${prefix}_prelim_coverage_filtered_diced.bed \
            > ${prefix}_max10kb_prelim_coverage_filtered_diced.bed
    
        """
    }

    process tfit_split_model_max5kb {
        println "[Log 4b]: Running Tfit model on <5kb prelim regions (maxk=2)"

        tag "$prefix"
        memory '70 GB'
        time '72h'
        queue 'long'
        clusterOptions = '-N 1 -n 32'

        publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*_bidir_predictions.bed"
        publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"

        when:
        params.tfit_split_model

        input:
        tuple val(prefix), file(bg), file(prelim_max5kb) from tfit_prelim_bg_out_split_max5kb

        output:
        tuple val(prefix), file ("*_bidir_predictions.bed") into tfit_model_bed_out_split_max5kb
        file ("*.tsv") into tfit_model_model_out_split_max5kb
        file ("*.log") into tfit_model_logs_out_split_max5kb

        script:
        """
        ${params.tfit_model_run} -t ${params.tfit_path} \
                            -c ${params.tfit_config_k2} \
                            -b ${bg} \
                            -k ${prelim_max5kb} \
                            -p ${prefix}_max5kb \
                            -n 32

        """
    }

    process tfit_split_model_max10kb {
        println "[Log 4b]: Running Tfit model on 5-10kb prelim regions (maxk=5)"

        tag "$prefix"
        memory '70 GB'
        time '96h'
        queue 'long'
        clusterOptions = '-N 1 -n 32'

        publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*_bidir_predictions.bed"
        publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"

        when:
        params.tfit_split_model

        input:
        tuple val(prefix), file(bg), file(prelim_max10kb) from tfit_prelim_bg_out_split_max10kb

        output:
        set val(prefix), file ("*_bidir_predictions.bed") into tfit_model_bed_out_split_max10kb
        file ("*.tsv") into tfit_model_model_out_split_max10kb
        file ("*.log") into tfit_model_logs_out_split_max10kb

        script:
        """
        ${params.tfit_model_run} -t ${params.tfit_path} \
                            -c ${params.tfit_config} \
                            -b ${bg} \
                            -k ${prelim_max10kb} \
                            -p ${prefix}_max10kb \
                            -n 32

        """
    }

    tfit_split_results = post_tfit_bg_split
        .join(tfit_model_bed_out_split_max5kb)
        .join(tfit_model_bed_out_split_max10kb)

    process tfit_split_cat {
        tag "$prefix"
        memory '4 GB'
        time '1h'
        queue 'short'

        publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*{split_bidir_predictions.bed,split_bidir_cov_filtered.bed}"

        when:
        params.tfit_split_model

        input:
        tuple val(prefix), file(bg), file(out_max5kb), file(out_max10kb) from tfit_split_results

        output:
        set val(prefix), file ("*split_bidir_predictions.bed"), file("*split_bidir_cov_filtered.bed") into tfit_model_out

        script:
        """
        grep -v '#' ${out_max5kb} > max5kb_noheader.bed
        grep -v '#' ${out_max10kb} > max10kb_noheader.bed
        cat max5kb_noheader.bed max10kb_noheader.bed | bedtools sort > ${prefix}_split_bidir_predictions.bed
        bedtools coverage -a ${prefix}_split_bidir_predictions.bed -b ${bg} > ${prefix}_split_bidir_cov.bed
        awk '{if (\$5 > 9) print \$0}' ${prefix}_split_bidir_cov.bed > ${prefix}_split_bidir_cov_filtered.bed
        """
    }

} else if (params.prelim_process && (params.tfit || params.tfit_model)) {

    tfit_prelim_bg_out_preprocess = prelimtfit_process_bg
        .join(tfit_prelim_out)

    process tfit_prelim_process {
        println "[Log 4b]: Processing prelim file"

        tag "$prefix"
        memory '20 GB'
        time '1h'
        queue 'short'

        publishDir "${params.outdir}/tfit/prelim", mode: 'copy', pattern: "*_prelim_coverage_filtered_diced.bed"

        when:
        params.prelim_process

        input:
        tuple val(prefix), file(bg), file(prelim) from tfit_prelim_bg_out_preprocess

        output:
        tuple val(prefix), file(bg), file ("*_prelim_coverage_filtered_diced.bed") into tfit_prelim_bg_out

        script:
        """
        python3 ${params.prelim_filter} \
                -p ${prelim} \
                -b ${bg} \
                -s ${prefix} \
                -o . \
                -g ${params.filtered_refseq} \
                -c ${params.chrom_sizes}

        """
    }

} else if (params.tfit || params.tfit_model) {

    tfit_prelim_bg_out = prelimtfit_process_bg
        .join(tfit_prelim_out)

}

if (params.tfit || params.tfit_model) {
    process tfit_model {
        println "[Log 4b]: Running Tfit model"

        tag "$prefix"
        memory '500 GB'
        time '72h'
        queue 'long'
        clusterOptions = '-N 1 -n 32'

        publishDir "${params.outdir}/tfit", mode: 'copy', pattern: "*{_bidir_predictions.bed,_bidir_cov_filtered.bed}"
        publishDir "${params.outdir}/tfit/logs", mode: 'copy', pattern: "*{tsv,log}"
    
        when:
        params.tfit_model || params.tfit

        input:
        tuple val(prefix), file(bg), file(prelim) from tfit_prelim_bg_out

        output:
        tuple val(prefix), file ("*_bidir_predictions.bed"), file("*_bidir_cov_filtered.bed") into tfit_model_out
        file ("*.tsv") into tfit_model_model_out
        file ("*.log") into tfit_model_logs_out

        script:
        """
        ${params.tfit_model_run} -t ${params.tfit_path} \
                            -c ${params.tfit_config} \
                            -b ${bg} \
                            -k ${prelim} \
                            -p ${prefix} \
                            -n 32

        bedtools coverage -a ${prefix}-1_bidir_predictions.bed -b ${bg} > ${prefix}_bidir_cov.bed
        awk '{if (\$5 > 9) print \$0}' ${prefix}_bidir_cov.bed > ${prefix}_bidir_cov_filtered.bed
        """
    }
}
println "[Log 4b]: Done Running Tfit model\n"


// PART 5: Preparing bigwig files for dREG

process dreg_prep {
    println "[Log 5]: Generating bigwig files for dREG"

    errorStrategy 'ignore'
    tag "$prefix"
    memory '60 GB'
    cpus 16
    queue 'short'

    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (params.savebg && (filename == "${prefix}.bedGraph"))       "bedgraphs/${prefix}.fivep.bedGraph"
             else if (params.savebw && (filename.indexOf("bw") > 0))        "bigwigs/$filename"
             else null
            }

    when:
    params.dreg || params.dreg_results

    input:
    set val(prefix), file(bam_file), file(index) from bam_for_dreg

    output:
    tuple val(prefix), file("${prefix}.pos.bw"), file("${prefix}.neg.bw") into dreg_bigwig
    tuple val(prefix), file("${prefix}.bedGraph") into dreg_bg

    script:
    if (params.singleEnd) {
        """
        echo "Creating BigWigs suitable as inputs to dREG"

        export CRAM_REFERENCE=${params.genome}

        bamToBed -i ${bam_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
        awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
        > ${prefix}.dreg.bed
        sortBed -i ${prefix}.dreg.bed > ${prefix}.dreg.sort.bed

        echo "positive strand processed to bedGraph"

        bedtools genomecov \
                -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${params.chrom_sizes} \
                -strand + \
                > ${prefix}.pos.bedGraph

        sortBed \
                -i ${prefix}.pos.bedGraph \
                > ${prefix}.pos.sort.bedGraph

        bedtools genomecov \
                -bg \
                -i ${prefix}.dreg.sort.bed \
                -g ${params.chrom_sizes} \
                -strand - \
                | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${prefix}.neg.bedGraph

        sortBed \
                -i ${prefix}.neg.bedGraph \
                > ${prefix}.neg.sort.bedGraph

        echo "negative strand processed to bedGraph"

        ${params.bedGraphToBigWig} ${prefix}.pos.sort.bedGraph ${params.chrom_sizes} ${prefix}.pos.bw
        ${params.bedGraphToBigWig} ${prefix}.neg.sort.bedGraph ${params.chrom_sizes} ${prefix}.neg.bw

        cat ${prefix}.pos.bedGraph \
        ${prefix}.neg.bedGraph \
        > ${prefix}.unsorted.bedGraph

        sortBed \
        -i ${prefix}.unsorted.bedGraph \
        > ${prefix}.bedGraph

        echo "bedGraph to bigwig done"
        """
    } else {
        if (params.r1_five_prime) {
            """
            samtools view -@ 16 -bf 0x2 ${bam_file} | samtools sort -n -@ 16 \
            > ${prefix}.dreg.bam

            bedtools bamtobed -bedpe -mate1 -i ${prefix}.dreg.bam \
              | awk 'BEGIN{OFS="\t"} (\$9 == "+") {print \$1,\$2,\$2+1,\$7,\$8,\$9}; (\$9 == "-") {print \$1,\$3-1,\$3,\$7,\$8,\$9}' \
              | sort -k 1,1 -k 2,2n > ${prefix}.dreg.sort.bed

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand + \
              > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand - \
              > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
              | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
              > ${prefix}.neg.bedGraph

            ${params.bedGraphToBigWig} ${prefix}.pos.bedGraph \
              ${params.chrom_sizes} ${prefix}.pos.bw

            ${params.bedGraphToBigWig} ${prefix}.neg.bedGraph \
              ${params.chrom_sizes} ${prefix}.neg.bw

            cat ${prefix}.pos.bedGraph \
            ${prefix}.neg.bedGraph \
            > ${prefix}.unsorted.bedGraph

            sortBed \
            -i ${prefix}.unsorted.bedGraph \
            > ${prefix}.bedGraph
            """
        } else {
            """
            samtools view -@ 16 -bf 0x2 ${bam_file} | samtools sort -n -@ 16 \
              > ${prefix}.dreg.bam

            bedtools bamtobed -bedpe -mate1 -i ${prefix}.dreg.bam \
              | awk 'BEGIN{OFS="\t"} (\$10 == "+") {print \$1,\$5,\$5+1,\$7,\$8,\$10}; (\$10 == "-") {print \$1,\$6-1,\$6,\$7,\$8,\$10}' \
              | sort -k 1,1 -k 2,2n > ${prefix}.dreg.sort.bed

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand + \
              > ${prefix}.pos.bedGraph

            bedtools genomecov -bg \
              -i ${prefix}.dreg.sort.bed \
              -g ${params.chrom_sizes} \
              -strand - \
              > ${prefix}.neg.noinv.bedGraph

            cat ${prefix}.neg.noinv.bedGraph \
              | awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,-1*\$4}' \
              > ${prefix}.neg.bedGraph

            ${params.bedGraphToBigWig} ${prefix}.pos.bedGraph \
              ${params.chrom_sizes} ${prefix}.pos.bw

            ${params.bedGraphToBigWig} ${prefix}.neg.bedGraph \
              ${params.chrom_sizes} ${prefix}.neg.bw

            cat ${prefix}.pos.bedGraph \
            ${prefix}.neg.bedGraph \
            > ${prefix}.unsorted.bedGraph

            sortBed \
            -i ${prefix}.unsorted.bedGraph \
            > ${prefix}.bedGraph
            """
        }
    }
}

println "[Log 5]: Bigwig files are ready \n"


// PART 6: Running dREG

process dreg_run {
    println "Log[6]: Running dREG"
    println "Log[6]: N.B. Requires GPUs"

    tag "$prefix"
    memory '50 GB'
    time '48h'
    cpus 4
    queue 'titan'
    clusterOptions '--gres=gpu'

    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "*dREG*"

    when:
    params.dreg

    input:
    tuple val(prefix), file(pos_bw), file(neg_bw) from dreg_bigwig

    output:
    tuple val(prefix), file ("${prefix}.*") into dREG_out

    script:
        """
        bash ${params.dreg_path} \
	     ${pos_bw} \
	     ${neg_bw} \
	     ${prefix} \
	     ${params.dreg_train} \
	     4 1 
        """
}

if (params.dreg_results) {
    dreg_res = Channel
        .fromPath(params.dreg_results)
        .map { file -> tuple((file.simpleName + '.sorted'), file)}

    dreg_res_process = dreg_bg
        .join(dreg_res)

} else if (params.dreg) {

    dreg_res_process = dreg_bg
        .join(dREG_out)

}

if (params.dreg_results || params.dreg) {
  process dreg_postprocess {
    println "Log[6]: Running dREG postprocessing"

    tag "$prefix"
    memory '8 GB'
    time '1h'
    cpus 1
    queue 'short'

    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "*covfiltered.bed"
    stageInMode 'copy'

    when:
    params.dreg || params.dreg_results

    input:
    tuple val(prefix), file(bg), file(dreg_resfile) from dreg_res_process

    output:
    tuple val(prefix), file ("${prefix}.dREG.full.covfiltered.bed") into dREG_res_out

    script:
        """
        gunzip ${prefix}.dREG.peak.full.bed.gz
	bedtools merge -i ${prefix}.dREG.peak.full.bed -d 20 > ${prefix}.dREG.peak.full.merge_distance20bp.bed
        bedtools coverage -a ${prefix}.dREG.peak.full.merge_distance20bp.bed -b ${bg} > ${prefix}.dREG.bidir.cov.bed
        awk '{if (\$4 > 9) print \$0}' ${prefix}.dREG.bidir.cov.bed > ${prefix}.dREG.full.covfiltered.bed
        gzip ${prefix}.dREG.peak.full.bed
        """
  }
}

// PART 7: Running stats on tfit/dreg outputs

if (params.savestats) {

  if (params.tfit || params.tfit_model || params.tfit_split_model) {

    println "[Log 7]: Getting tfit bidir stats"

    process calculate_tfit_stats {
      tag "$prefix"
      memory '8 GB'
      time '1h'
      cpus 1
      queue 'short'
      stageInMode 'copy'

      publishDir "${params.outdir}" , mode: 'copy',
      saveAs: {filename ->
               if (params.savebidirs && (filename.indexOf(".tfit.") > 0))       "tfit_subsets/$filename"
               else null
              }

      when:
      params.tfit || params.tfit_split_model || params.tfit_model

      input:
      tuple val(prefix), file(nonfilt), file(tfit_regions) from tfit_model_out

      output:
      tuple val(prefix), file ("${prefix}.tfit.*") into tfit_subsets
      file ("${prefix}.tfit_stats.txt") into tfit_stats

      script:
      """
      bedtools intersect -wa -a ${tfit_regions} -b ${params.promoter_bed} > ${prefix}.tfit.promoter_regions.bed
      bedtools intersect -wa -a ${tfit_regions} -b ${params.intron_bed} > intron_regions_temp.bed
      bedtools intersect -wa -a ${tfit_regions} -b ${params.intergenic_bed} > intergenic_regions_temp.bed
      bedtools intersect -v -a intron_regions_temp.bed -b ${params.exon_bed} \
        > intron_regions_temp2.bed
      bedtools intersect -v -a intron_regions_temp2.bed -b ${prefix}.tfit.promoter_regions.bed \
        > ${prefix}.tfit.intron_regions.bed
      bedtools intersect -v -a intergenic_regions_temp.bed -b ${params.exon_bed} \
        > intergenic_regions_temp2.bed
      bedtools intersect -v -a intergenic_regions_temp2.bed -b ${prefix}.tfit.promoter_regions.bed \
        > ${prefix}.tfit.intergenic_regions.bed

      bedtools getfasta -fi ${params.genome} -bed ${tfit_regions} -fo seqs.txt
      gcnum=\$(tr -cd 'CGcg' < seqs.txt | wc -c)
      totnum=\$(tr -cd 'ATCGatcg' < seqs.txt | wc -c)
      gc_prop=\$(echo "scale=3 ; \$gcnum / \$totnum" | bc)

      printf "%s\t%s\t%s\t%s\t%s\t%.3f\n" \
          ${prefix} \
          \$(wc -l ${tfit_regions} | awk '{print \$1}') \
          \$(wc -l ${prefix}.tfit.promoter_regions.bed | awk '{print \$1}') \
          \$(wc -l ${prefix}.tfit.intron_regions.bed | awk '{print \$1}') \
          \$(wc -l ${prefix}.tfit.intergenic_regions.bed | awk '{print \$1}') \
          \$gc_prop \
          > ${prefix}.tfit_stats.txt
      """
    }

    process accumulate_tfit_stats {
      tag "tfit_stats"
      memory '1 GB'
      time '1h'
      cpus 1
      queue 'short'
      stageInMode 'copy'

      publishDir "${params.outdir}/bidir_summary" , mode: 'copy', pattern: "tfit_stats.txt"

      when:
      params.tfit || params.tfit_split_model || params.tfit_model

      input:
      file ('*') from tfit_stats.collect()

      output:
      file ("tfit_stats.txt") into tfit_stats_accum

      script:
      """
      printf "sample_name\tnum_tfit_bidir\tnum_tfit_bidir_promoter\tnum_tfit_bidir_intronic\tnum_tfit_bidir_intergenic\ttfit_bidir_gc_prop\n" \
      > tfit_stats.txt
      cat *.tfit_stats.txt >> tfit_stats.txt

      """
    }

  } else if (params.dreg || params.dreg_results) {

    process calculate_dreg_stats {
      tag "$prefix"
      memory '8 GB'
      time '1h'
      cpus 1
      queue 'short'
      stageInMode 'copy'

      publishDir "${params.outdir}" , mode: 'copy',
      saveAs: {filename ->
               if (params.savebidirs && (filename.indexOf(".dreg.") > 0))       "dreg_subsets/$filename"
               else null
              }

      when:
      params.dreg || params.dreg_results

      input:
      tuple val(prefix), file(dreg_regions) from dREG_res_out

      output:
      tuple val(prefix), file ("${prefix}.dreg.*") into dreg_subsets
      file ("${prefix}.dreg_stats.txt") into dreg_stats

      script:
      """
      bedtools intersect -wa -a ${dreg_regions} -b ${params.promoter_bed} > ${prefix}.dreg.promoter_regions.bed
      bedtools intersect -wa -a ${dreg_regions} -b ${params.intron_bed} > intron_regions_temp.bed
      bedtools intersect -wa -a ${dreg_regions} -b ${params.intergenic_bed} > intergenic_regions_temp.bed
      bedtools intersect -v -a intron_regions_temp.bed -b ${params.exon_bed} \
        > intron_regions_temp2.bed
      bedtools intersect -v -a intron_regions_temp2.bed -b ${prefix}.dreg.promoter_regions.bed \
        > ${prefix}.dreg.intron_regions.bed
      bedtools intersect -v -a intergenic_regions_temp.bed -b ${params.exon_bed} \
        > intergenic_regions_temp2.bed
      bedtools intersect -v -a intergenic_regions_temp2.bed -b ${prefix}.dreg.promoter_regions.bed \
        > ${prefix}.dreg.intergenic_regions.bed

      bedtools getfasta -fi ${params.genome} -bed ${dreg_regions} -fo seqs.txt
      gcnum=\$(tr -cd 'CGcg' < seqs.txt | wc -c)
      totnum=\$(tr -cd 'ATCGatcg' < seqs.txt | wc -c)
      gc_prop=\$(echo "scale=3 ; \$gcnum / \$totnum" | bc)

      printf "%s\t%s\t%s\t%s\t%s\t%.3f\n" \
          ${prefix} \
          \$(wc -l ${dreg_regions} | awk '{print \$1}') \
          \$(wc -l ${prefix}.dreg.promoter_regions.bed | awk '{print \$1}') \
          \$(wc -l ${prefix}.dreg.intron_regions.bed | awk '{print \$1}') \
          \$(wc -l ${prefix}.dreg.intergenic_regions.bed | awk '{print \$1}') \
          \$gc_prop \
          > ${prefix}.dreg_stats.txt
      """
    }

    process accumulate_dreg_stats {
      tag "dreg_stats"
      memory '1 GB'
      time '1h'
      cpus 1
      queue 'short'
      stageInMode 'copy'

      publishDir "${params.outdir}/bidir_summary" , mode: 'copy', pattern: "dreg_stats.txt"

      when:
      params.dreg || params.dreg_results

      input:
      file ('*') from dreg_stats.collect()

      output:
      file ("dreg_stats.txt") into dreg_stats_accum

      script:
      """
      printf "sample_name\tnum_tfit_bidir\tnum_tfit_bidir_promoter\tnum_tfit_bidir_intronic\tnum_tfit_bidir_intergenic\ttfit_bidir_gc_prop\n" \
      > dreg_stats.txt
      cat *.dreg_stats.txt >> dreg_stats.txt

      """
    }
    println "[Log 7]: Finished with bidir stats"

  } else {
    null

  }
}

// PART 8: Counting over genes

process gene_count {
   println "[Log 8]: Running FeatureCounts"

    tag "$prefix"
    memory '8 GB'
    time '3h'
    cpus 8
    queue 'short'

    publishDir "${params.outdir}/featurecounts_genes/", mode: 'copy', pattern: "*gene_counts.txt"

    when:
    params.gene_count

    input:
    tuple val(prefix), file(bam_file), file(index) from bam_for_gene_counting

    output:
    tuple val(prefix), file ("*gene_counts.txt") into gene_count_out

    script:
    if (params.singleEnd) {
        paired = 'FALSE'
    } else {
        paired = 'TRUE'
    }
    
    """
    #!/usr/bin/env Rscript

    library("Rsubread")

    gtf_table <- read.table("${params.filtered_refseq}")   

    if (${paired} == 'FALSE') {

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=1,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.trunc_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=1,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.5ptrunc_gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    } else {

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=2,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.trunc_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=2,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".stranded.5ptrunc_gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    }

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.filtered_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=0,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".unstranded.gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE) 

    gtf_table <- read.table("${params.trunc_refseq}")
    
    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.trunc_refseq}",
        isGTFAnnotationFile=TRUE,
        GTF.featureType="gene_length",
        useMetaFeatures=FALSE,
        allowMultiOverlap=TRUE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=0,
        nthreads=8)
    fc\$annotation["TranscriptID"] <- gtf_table["V13"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","TranscriptID","Length")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".unstranded.5ptrunc_gene_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)
    
    """
}

println "[Log 8]: Done Running FeatureCounts\n"

// PART 9: Counting over bidirectionals

process bidirectional_count {
    println "[Log 9]: Running FeatureCounts over bidirs"

    tag "$prefix"
    memory '8 GB'
    time '3h'
    cpus 8
    queue 'short'

    publishDir "${params.outdir}/featurecounts_bidirs/", mode: 'copy', pattern: "*bidir_counts.txt"

    when:
    params.bidir_count

    input:
    tuple val(prefix), file(bam_file), file(index) from bam_for_bidir_counting

    output:
    tuple val(prefix), file ("*bidir_counts.txt") into bidir_count_out

    script:
    if (params.singleEnd) {
        paired = 'FALSE'
    } else {
        paired = 'TRUE'
    }

    """
    #!/usr/bin/env Rscript

    library("Rsubread")

    saf_table <- read.table("${params.bidir_accum}", header=TRUE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.bidir_accum}",
        isGTFAnnotationFile=FALSE,
        useMetaFeatures=FALSE,
        allowMultiOverlap=FALSE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=1,
        nthreads=8)
    fc\$annotation["Source"] <- saf_table["Source"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","Source")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".pos.bidir_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.bidir_accum}",
        isGTFAnnotationFile=FALSE,
        useMetaFeatures=FALSE,
        allowMultiOverlap=FALSE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=2,
        nthreads=8)
    fc\$annotation["Source"] <- saf_table["Source"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","Source")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".neg.bidir_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    fc <- featureCounts(files="${bam_file}",
        annot.ext="${params.bidir_accum}",
        isGTFAnnotationFile=FALSE,
        useMetaFeatures=FALSE,
        allowMultiOverlap=FALSE,
        largestOverlap=TRUE,
        countMultiMappingReads=FALSE,
        isPairedEnd=${paired},
        strandSpecific=0,
        nthreads=8)
    fc\$annotation["Source"] <- saf_table["Source"]
    write.table(x=data.frame(fc\$annotation[,c("GeneID","Source")],
                             fc\$counts,stringsAsFactors=FALSE),
        file=paste0("${prefix}",".unstranded.bidir_counts.txt"),
        quote=FALSE,sep="\t",
        row.names=FALSE)

    """

}

println "[Log 9]: Done Running FeatureCounts\n"

/*
 * Completion report
 */
workflow.onComplete {

    def report_fields = [:]
    report_fields['version'] = params.version
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    report_fields['summary']['Pipeline repository Git URL'] = workflow.repository ?: 'Not stored'
    report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId ?: 'See hash'
    report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision ?: 'See hash'
    report_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    report_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    report_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/report_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/report_template.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report_bidir_${workflow.runName}.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report_bidir_${workflow.runName}.txt" )
    output_tf.withWriter { w -> w << report_txt }

    log.info "[Bidirectional-Flow] Pipeline Complete"

}

