include { BCFTOOLS_NORM   } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_VIEW   } from '../../../modules/nf-core/bcftools/view'
include { VCFLIB_VCFFIXUP } from '../../../modules/nf-core/vcflib/vcffixup/main'

workflow VCF_NORMALIZE_BCFTOOLS {
    take:
    ch_vcf_tbi      // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, [fai, gzi] ]
    normalize       // boolean
    compute_freq    // boolean

    main:

    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Join duplicated biallelic sites into multiallelic records
    if (normalize) {
        BCFTOOLS_NORM(ch_vcf_tbi, ch_fasta)

        // Join multiallelic VCF and TBI
        ch_multiallelic_vcf_tbi = BCFTOOLS_NORM.out.vcf
            .join(
                BCFTOOLS_NORM.out.tbi.mix(
                    BCFTOOLS_NORM.out.csi
                )
            )

        // Remove all multiallelic records and samples specified in the `--remove_samples` command:
        BCFTOOLS_VIEW(ch_multiallelic_vcf_tbi, [], [], [])

        // Join biallelic VCF and TBI
        ch_vcf_tbi = BCFTOOLS_VIEW.out.vcf
            .join(
                BCFTOOLS_VIEW.out.tbi.mix(
                    BCFTOOLS_VIEW.out.csi
                )
            )
    }

    // (Optional) Fix panel (When AC/AN INFO fields in VCF are inconsistent with GT field)
    if (compute_freq) {
        VCFLIB_VCFFIXUP(ch_vcf_tbi)

        // Join fixed vcf and tbi
        ch_vcf_tbi = VCFLIB_VCFFIXUP.out.vcf
            .join(VCFLIB_VCFFIXUP.out.index)
    }
    emit:
    vcf_tbi        = ch_vcf_tbi                     // channel: [ [id, chr], vcf, tbi ]
}
