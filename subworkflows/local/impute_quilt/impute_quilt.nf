include { QUILT_QUILT              } from '../../../modules/nf-core/quilt/quilt'
include { BCFTOOLS_ANNOTATE        } from '../../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'


workflow BAM_IMPUTE_QUILT {

    take:
    ch_hap_legend        // channel: [ [panel, chr], hap, legend ]
    ch_input             // channel: [ [id, chr], bam, bai ]
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
        ch_hap_chunks = ch_hap_legend.join(ch_chunks, by:0).join(genetic_map_file)
    }

    ch_quilt = ch_input
        .map{ metaIC, bam, bai -> [metaIC.subMap("chr"), metaIC, bam, bai]}
        .combine(ch_hap_chunks
            .map{ metaIC, hap, legend, chr, start, end, ngen, buffer, gmap ->
                [metaIC.subMap("chr"), metaIC, hap, legend, chr, start, end, ngen, buffer, gmap]
            }, by:0
        )
        .map {
            metaC, metaIC, bam, bai, metaPC, hap, legend, chr, start, end, ngen, buffer, gmap ->
            [metaIC + ["panel": metaPC.id], bam, bai, hap, legend, chr, start, end, ngen, buffer, gmap]
        }

    // Run QUILT
    QUILT_QUILT ( ch_quilt, posfile_phasefile, fasta )
    ch_versions = ch_versions.mix(QUILT_QUILT.out.versions.first())

    // Index imputed VCF
    BCFTOOLS_INDEX_1(QUILT_QUILT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())

    // Annotate the variants
    BCFTOOLS_ANNOTATE(QUILT_QUILT.out.vcf
        .join(BCFTOOLS_INDEX_1.out.tbi)
        .combine(Channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    // Index imputed annotated VCF
    BCFTOOLS_INDEX_2(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())

    // Join VCFs and TBIs
    ch_vcf_tbi = BCFTOOLS_ANNOTATE.out.vcf.join(BCFTOOLS_INDEX_2.out.tbi)

    emit:
    vcf_tbi     = ch_vcf_tbi               // channel:  [ meta, vcf, tbi ]
    versions    = ch_versions              // channel:   [ versions.yml ]
}
