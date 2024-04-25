include { STITCH             } from '../../../modules/nf-core/stitch/main'
include { BCFTOOLS_INDEX     } from '../../../modules/nf-core/bcftools/index/main'


workflow BAM_IMPUTE_STITCH {

    take:
    ch_parameters                                //  channel: [ val(meta), bam, bai ]
    ch_samples
    ch_fasta

    main:

    ch_versions      = Channel.empty()

    // Run STITCH
    seed = params.seed
    STITCH( ch_parameters, ch_samples, ch_fasta, seed )

    // Index imputed annotated VCF
    BCFTOOLS_INDEX(STITCH.out.vcf)

    // Join VCFs and TBIs
    ch_vcf_tbi = STITCH.out.vcf.join(BCFTOOLS_INDEX.out.tbi)


    emit:
    vcf_tbi                    = ch_vcf_tbi                            // channel:   [ meta, vcf, tbi ]
    versions                   = ch_versions                           // channel:   [ versions.yml ]

}
