include { STITCH               } from '../../../modules/nf-core/stitch'
include { BCFTOOLS_INDEX       } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_STITCH {

    take:
    ch_parameters  // channel:   [ [chr], posfile, input, rdata, chr, k_val, ngen]
    ch_bam         // channel:   [ [id], bam, bai, bamlist ]
    ch_fasta       // channel:   [ [genome], fa, fai ]

    main:

    ch_versions      = Channel.empty()
    // Run STITCH
    seed = params.seed

    ch_bam_params = ch_bam // Add chr to meta map
        .combine(ch_parameters)
        .map{
            metaI, bam, bai, bamlist, metaPC, posfile, input, rdata, chr, k_val, ngen ->
            [metaI + [chr: metaPC.chr, panel:metaPC.id], bam, bai, bamlist, posfile, input, rdata, chr, k_val, ngen]
        }

    STITCH( ch_bam_params, ch_fasta, seed )
    ch_versions = ch_versions.mix(STITCH.out.versions)

    // Index imputed annotated VCF
    BCFTOOLS_INDEX(STITCH.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    // Join VCFs and TBIs
    ch_vcf_tbi = STITCH.out.vcf
        .join(BCFTOOLS_INDEX.out.tbi)
        .map { metaI, vcf, tbi -> [ metaI + [tools: "stitch"], vcf, tbi ] }

    emit:
    vcf_tbi  = ch_vcf_tbi                        // channel:   [ [id, chr], vcf, tbi ]
    versions = ch_versions                       // channel:   [ versions.yml ]

}
