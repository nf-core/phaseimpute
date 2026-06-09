include { BCFTOOLS_NORM   } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_VIEW   } from '../../../modules/nf-core/bcftools/view'
include { VCFLIB_VCFFIXUP } from '../../../modules/nf-core/vcflib/vcffixup/main'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index'

workflow VCF_NORMALIZE_BCFTOOLS {
    take:
    ch_vcf_index    // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, [fai, gzi] ]
    normalize       // boolean
    compute_freq    // boolean

    main:

    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Join duplicated biallelic sites into multiallelic records
    if (normalize) {
        BCFTOOLS_NORM(ch_vcf_index, ch_fasta)

        // Join multiallelic VCF and TBI
        ch_multiallelic_vcf_index = BCFTOOLS_NORM.out.vcf
            .join(BCFTOOLS_NORM.out.index)

        // Remove all multiallelic records and samples specified in the `--remove_samples` command:
        BCFTOOLS_VIEW(ch_multiallelic_vcf_index, [], [], [])

        // Join biallelic VCF and index
        ch_vcf_index = BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_VIEW.out.index)
    }

    // (Optional) Fix panel (When AC/AN INFO fields in VCF are inconsistent with GT field)
    if (compute_freq) {
        VCFLIB_VCFFIXUP(ch_vcf_index)

        // Index fixed panel
        BCFTOOLS_INDEX(VCFLIB_VCFFIXUP.out.vcf)

        // Join fixed VCF and index
        ch_vcf_index = VCFLIB_VCFFIXUP.out.vcf
            .join(BCFTOOLS_INDEX.out.index)
    }
    emit:
    vcf_index = ch_vcf_index // channel: [ [id, chr], vcf, [tbi, csi] ]
}
