include { QUILT_QUILT                        } from '../../../modules/nf-core/quilt/quilt'
include { BCFTOOLS_ANNOTATE                  } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_IMPUTE_QUILT {

    take:
    ch_input             // channel: [ [id], bam, bai ]
    ch_hap_legend        // channel: [ [panel, chr], hap, legend ]
    ch_chunks            // channel: [ [panel, chr], chr, start_coordinate, end_coordinate ]
    ch_fasta             // channel: [ [genome], fa, fai ]

    main:

    ch_versions = Channel.empty()

    posfile             = []
    phasefile           = []
    posfile_phasefile   = [[id: null], posfile, phasefile]
    genetic_map_file    = []

    ngen                = params.ngen
    buffer              = params.buffer

    ch_hap_chunks = ch_hap_legend
        .combine(ch_chunks, by:0)
        .map { it + ngen + buffer + [[]] }

    if (!genetic_map_file.isEmpty()) {
        // Add genetic map file (untested)
        ch_hap_chunks = ch_hap_chunks
            .map{it[0..-1]}
            .join(genetic_map_file)
    }

    ch_quilt = ch_input
        .map{ metaI, bam, bai -> [[id: "all"], metaI, bam, bai] }
        .groupTuple()
        .map { metaI, all_metas, bam, bai -> [metaI + [metas: all_metas], bam, bai] }
        .combine(ch_hap_chunks)
        .map {
            metaI, bam, bai, metaPC, hap, legend, chr, start, end, ngen, buffer, gmap ->
            [
                metaI + [panel: metaPC.id, chr: metaPC.chr, chunk: start + "-" + end],
                bam, bai, hap, legend, chr, start, end, ngen, buffer, gmap
            ]
        }

    // Run QUILT
    QUILT_QUILT ( ch_quilt, posfile_phasefile, ch_fasta )
    ch_versions = ch_versions.mix(QUILT_QUILT.out.versions.first())

    // Annotate the variants
    BCFTOOLS_ANNOTATE(QUILT_QUILT.out.vcf
        .join(QUILT_QUILT.out.tbi)
        .combine(Channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    // Join VCFs and TBIs
    ch_vcf_tbi = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi)
        .map { metaIPC, vcf, tbi -> [metaIPC + [tools: "quilt"], vcf, tbi] }

    emit:
    vcf_tbi     = ch_vcf_tbi               // channel:  [ [id, panel], vcf, tbi ]
    versions    = ch_versions              // channel:  [ versions.yml ]
}
