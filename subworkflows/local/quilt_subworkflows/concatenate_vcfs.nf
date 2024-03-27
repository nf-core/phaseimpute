include { BCFTOOLS_INDEX  } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'

workflow CONCATENATE_VCF {

    take:
    ch_imputedvcf                            // channel: [ val(meta),bam ]

    main:

    // Index imputed VCF
    BCFTOOLS_INDEX(ch_imputedvcf)

    // Join VCFs and TBIs
    ch_vcf_tbi = ch_imputedvcf.join(BCFTOOLS_INDEX.out.tbi)
    ch_vcf_tbi.dump(tag:"ch_vcf_tbi")

    // Remove chromosome from meta
    ch_vcf_tbi_grouped = ch_vcf_tbi.map{ meta, vcf, tbi ->
                        return [['id' : meta.id], vcf, tbi]
                        }

    ch_vcf_tbi_grouped = ch_vcf_tbi_grouped.groupTuple(
        by:[0],
        sort:
            { a, b ->
                def fa = a.name.tokenize('.')
                def fb = b.name.tokenize('.')

                // Compare chr  ?: Compare start ?: Compare stop
                fa[1] <=> fb[1] ?: fa[2] <=> fb[2] ?: fa[3] <=> fb[3]
            }
    )


    // Ligate and concatenate chunks
    BCFTOOLS_CONCAT(ch_vcf_tbi_grouped)

    emit:
    ch_concat_vcf =  BCFTOOLS_CONCAT.out.vcf

    }
