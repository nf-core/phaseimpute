include { QUILT_QUILT                        } from '../../../modules/nf-core/quilt/quilt'
include { BCFTOOLS_ANNOTATE                  } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_IMPUTE_QUILT {

    take:
    ch_input             // channel: [ [id], [bam], [bai], bampaths, bamnames ]
    ch_hap_legend        // channel: [ [panel_id, chr], hap, legend ]
    ch_chunks            // channel: [ [panel_id, chr], chr, start_coordinate, end_coordinate ]
    ch_fasta             // channel: [ [genome], fa, fai ]

    main:

    ch_versions = Channel.empty()

    genetic_map_file    = []

    ngen_params         = params.ngen
    buffer_params       = params.buffer

    ch_hap_chunks = ch_hap_legend
        .combine(ch_chunks, by:0)
        .map {it -> it + ngen_params + buffer_params + [[]] }

    if (!genetic_map_file.isEmpty()) {
        // Add genetic map file (untested)
        ch_hap_chunks = ch_hap_chunks
            .map{it -> it[0..-1]}
            .join(genetic_map_file)
    }

    ch_quilt = ch_input
        .combine(ch_hap_chunks)
        .map {
            metaI, bam, bai, bampath, bamnames, metaPC, hap, legend, chr, start, end, ngen, buffer, gmap ->
            [
                metaI + [panel_id: metaPC.panel_id, chr: metaPC.chr, chunk: chr + ":" + start + "-" + end],
                bam, bai, bampath, bamnames, hap, legend, [], [], [], chr, start, end, ngen, buffer, gmap
            ]
        }

    // Run QUILT
    QUILT_QUILT ( ch_quilt, ch_fasta )
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
    vcf_tbi     = ch_vcf_tbi               // channel:  [ [id, panel_id, tools], vcf, tbi ]
    versions    = ch_versions              // channel:  [ versions.yml ]
}
