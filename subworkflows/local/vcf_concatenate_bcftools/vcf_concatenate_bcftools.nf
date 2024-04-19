include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'

workflow VCF_CONCATENATE_BCFTOOLS {

    take:
    ch_vcf_tbi                            // channel: [ val(meta), vcf, tbi ]

    main:

    // Remove chromosome from meta
    ch_vcf_tbi_grouped = ch_vcf_tbi.map{ meta, vcf, tbi ->
                        return [['id' : meta.id], vcf, tbi]
                        }
    // Group by ID
    ch_vcf_tbi_grouped = ch_vcf_tbi_grouped.groupTuple( by:[0] )

    // Ligate and concatenate chunks
    BCFTOOLS_CONCAT(ch_vcf_tbi_grouped)

    emit:
    ch_concat_vcf =  BCFTOOLS_CONCAT.out.vcf

    }
