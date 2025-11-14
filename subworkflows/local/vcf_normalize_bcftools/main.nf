include { BCFTOOLS_NORM   } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_VIEW   } from '../../../modules/nf-core/bcftools/view'
include { VCFLIB_VCFFIXUP } from '../../../modules/nf-core/vcflib/vcffixup/main'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index'

workflow VCF_NORMALIZE_BCFTOOLS {
    take:
    ch_vcf_tbi      // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()
    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Join duplicated biallelic sites into multiallelic records
    if (params.normalize) {
        BCFTOOLS_NORM(ch_vcf_tbi, ch_fasta)
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

        // Join multiallelic VCF and TBI
        ch_multiallelic_vcf_tbi = BCFTOOLS_NORM.out.vcf
            .join(BCFTOOLS_NORM.out.tbi)

        // Remove all multiallelic records and samples specified in the `--remove_samples` command:
        BCFTOOLS_VIEW(ch_multiallelic_vcf_tbi, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

        // Join biallelic VCF and TBI
        ch_vcf_tbi = BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_VIEW.out.tbi)
    }

    // (Optional) Fix panel (When AC/AN INFO fields in VCF are inconsistent with GT field)
    if (params.compute_freq == true) {
        VCFLIB_VCFFIXUP(ch_vcf_tbi)
        ch_versions = ch_versions.mix(VCFLIB_VCFFIXUP.out.versions.first())

        // Index fixed panel
        BCFTOOLS_INDEX(VCFLIB_VCFFIXUP.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

        // Join fixed vcf and tbi
        ch_vcf_tbi = VCFLIB_VCFFIXUP.out.vcf
            .join(BCFTOOLS_INDEX.out.tbi)
    }
    emit:
    vcf_tbi        = ch_vcf_tbi                     // channel: [ [id, chr], vcf, tbi ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}
