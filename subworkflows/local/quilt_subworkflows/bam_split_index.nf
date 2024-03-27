include { BAMTOOLS_SPLIT } from '../../modules/nf-core/bamtools/split/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow BAM_SPLIT_INDEX {

    take:
    ch_single_bam                             // channel: [ val(meta),bam ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: BAMTOOLS_SPLIT
    //

    BAMTOOLS_SPLIT( ch_single_bam )
    ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions)

    ch_transpose_bams = BAMTOOLS_SPLIT.out.bam.transpose()

    ch_numeric_bams = ch_transpose_bams.flatMap { meta, file ->
        def matcher = file =~ /.*\.REF_(\d+|X)\.bam$/
        if (matcher.matches()) {
            def chr = matcher[0][1]
            return [[['id':meta.id,'chr': chr], file]]
        } else {
            return []
            }
        }

    //
    // MODULE: SAMTOOLS_INDEX
    //

    SAMTOOLS_INDEX( ch_numeric_bams )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_numeric_bams_bai = ch_numeric_bams.join(SAMTOOLS_INDEX.out.bai)

    ch_numeric_bams_bai = ch_numeric_bams_bai.map{ meta, bam, bai ->
                    [meta.chr, meta, bam, bai]
                    }


    emit:
    ch_numeric_bams_bai       = ch_numeric_bams_bai                    // channel: [ chr, val(meta), bam, bai ]
    versions                  = ch_versions                           // channel: [ versions.yml ]
}
