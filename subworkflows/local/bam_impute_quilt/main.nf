include { QUILT_QUILT                        } from '../../../modules/nf-core/quilt/quilt'
include { BCFTOOLS_INDEX                     } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_ANNOTATE                  } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_IMPUTE_QUILT {

    take:
    ch_input             // channel: [ [id], bam, bai ]
    ch_hap_legend        // channel: [ [panel, chr], hap, legend ]
    ch_chunks            // channel: [ [panel, chr], start_coordinate, end_coordinate, number ]


    main:


    ch_versions = Channel.empty()

    posfile             = []
    phasefile           = []
    posfile_phasefile   = [[id: null], posfile, phasefile]
    genetic_map_file    = []
    fasta               = [[id:'test'], []]

    ngen                = params.ngen
    buffer              = params.buffer

    if (genetic_map_file.isEmpty()) {
        ch_hap_chunks = ch_hap_legend.combine(ch_chunks, by:0).map { it + ngen + buffer + [[]] }
    } else {
        // Add ngen and buffer + genetic map file (untested)
        ch_hap_chunks = ch_hap_legend.combine(ch_chunks, by:0).join(genetic_map_file)
    }

    ch_quilt = ch_input
        .combine(ch_hap_chunks)
        .map {
            metaIC, bam, bai, metaPC, hap, legend, chr, start, end, ngen, buffer, gmap ->
            [
                metaIC.subMap("id") + ["panel": metaPC.id, "chr": metaPC.chr, "chunk": start + "-" + end],
                bam, bai, hap, legend, chr, start, end, ngen, buffer, gmap
            ]
        }

    // Run QUILT
    QUILT_QUILT ( ch_quilt, posfile_phasefile, fasta )
    ch_versions = ch_versions.mix(QUILT_QUILT.out.versions.first())

    // Index imputed VCF
    BCFTOOLS_INDEX(QUILT_QUILT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    // Annotate the variants
    BCFTOOLS_ANNOTATE(QUILT_QUILT.out.vcf
        .join(BCFTOOLS_INDEX.out.tbi)
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
