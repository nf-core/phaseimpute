include { QUILT_QUILT     } from '../../../modules/nf-core/quilt/quilt/main'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index/main'


workflow IMPUTE_QUILT {

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

    ch_quilt = ch_input.combine(ch_hap_chunks)
    ch_quilt_input = ch_quilt.map { it.take(3) + it.drop(4) }

    // Add metamap with chromosome information
    ch_quilt_input = ch_quilt_input
        .map{ meta, bam, bai, hap, legend, chr, start, end, ngen2, buffer2, genetic ->
            [['id': meta.id, 'chr': chr], bam, bai, hap, legend, chr, start, end, ngen2, buffer2, genetic]
        }

    // Run QUILT
    QUILT_QUILT ( ch_quilt_input, posfile_phasefile, fasta )

    // Index imputed VCF
    BCFTOOLS_INDEX(QUILT_QUILT.out.vcf)

    // Join VCFs and TBIs
    ch_vcf_tbi = QUILT_QUILT.out.vcf.join(BCFTOOLS_INDEX.out.tbi)

    emit:
    vcf_tbi     = ch_vcf_tbi               // channel:  [ meta, vcf, tbi ]
    versions    = ch_versions              // channel:   [ versions.yml ]
}
