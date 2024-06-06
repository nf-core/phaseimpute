include { BCFTOOLS_NORM                         } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_3    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_4    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_VIEW as BCFTOOLS_DEL_MLT_ALL } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_DEL_SPL     } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_CONVERT                      } from '../../../modules/nf-core/bcftools/convert'
include { VCFLIB_VCFFIXUP                       } from '../../../modules/nf-core/vcflib/vcffixup/main'


workflow VCF_NORMALIZE_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()
    ch_fasta = ch_fasta.map { meta, fasta, fai -> [meta, fasta] }

    // Join duplicated biallelic sites into multiallelic records
    BCFTOOLS_NORM(ch_vcf, ch_fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // Index multiallelic VCF
    BCFTOOLS_INDEX_1(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions)

    // Join multiallelic VCF and TBI
    ch_multiallelic_vcf_tbi = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_INDEX_1.out.tbi)

    // Remove all multiallelic records:
    BCFTOOLS_DEL_MLT_ALL(ch_multiallelic_vcf_tbi, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_DEL_MLT_ALL.out.versions)

    // Index biallelic VCF
    BCFTOOLS_INDEX_2(BCFTOOLS_DEL_MLT_ALL.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions)

    // Join biallelic VCF and TBI
    ch_biallelic_vcf_tbi = BCFTOOLS_DEL_MLT_ALL.out.vcf.join(BCFTOOLS_INDEX_2.out.tbi)

    // (Optional) Remove benchmarking samples (e.g. NA12878) from the reference panel
    if (!(params.remove_samples == null)){
        BCFTOOLS_DEL_SPL(ch_biallelic_vcf_tbi, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_DEL_SPL.out.versions)

        BCFTOOLS_INDEX_3(BCFTOOLS_DEL_SPL.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_3.out.versions)

        ch_biallelic_vcf_tbi_spl = BCFTOOLS_DEL_SPL.out.vcf.join(BCFTOOLS_INDEX_3.out.tbi)
    } else {
        ch_biallelic_vcf_tbi_spl = ch_biallelic_vcf_tbi
    }

    if (params.compute_freq == true) {
        // Fix panel (AC/AN INFO fields in VCF are inconsistent with GT field)
        VCFLIB_VCFFIXUP(ch_biallelic_vcf_tbi_spl)
        ch_versions = ch_versions.mix(VCFLIB_VCFFIXUP.out.versions)

        // Index fixed panel
        BCFTOOLS_INDEX_4(VCFLIB_VCFFIXUP.out.vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX_4.out.versions)

        // Join fixed vcf and tbi
        ch_biallelic_vcf_tbi_freq = VCFLIB_VCFFIXUP.out.vcf.join(BCFTOOLS_INDEX_4.out.tbi)
    } else {
        ch_biallelic_vcf_tbi_freq = ch_biallelic_vcf_tbi_spl
    }

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_biallelic_vcf_tbi_freq, ch_fasta, [])
    ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions)

    // Output hap and legend files
    ch_hap_legend = BCFTOOLS_CONVERT.out.hap.join(BCFTOOLS_CONVERT.out.legend)

    emit:
    vcf_tbi        = ch_biallelic_vcf_tbi_freq      // channel: [ [id, chr], vcf, tbi ]
    hap_legend     = ch_hap_legend                  // channel: [ [id, chr], '.hap', '.legend' ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}
