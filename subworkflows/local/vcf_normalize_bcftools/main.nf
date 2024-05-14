include { BCFTOOLS_NORM                         } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_3    } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_VIEW as BCFTOOLS_DEL_MLT_ALL } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_DEL_SPL     } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_CONVERT                      } from '../../../modules/nf-core/bcftools/convert'


workflow VCF_NORMALIZE_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()
    ch_fasta = ch_fasta.map { meta, fasta, fai -> [meta, fasta] }

    // Join duplicated biallelic sites into multiallelic records
    BCFTOOLS_NORM(ch_vcf, ch_fasta)

    // Index multiallelic VCF
    BCFTOOLS_INDEX_1(BCFTOOLS_NORM.out.vcf)

    // Join multiallelic VCF and TBI
    ch_multiallelic_vcf_tbi = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_INDEX_1.out.tbi)

    // Remove all multiallelic records:
    BCFTOOLS_DEL_MLT_ALL(ch_multiallelic_vcf_tbi, [], [], [])

    // Index biallelic VCF
    BCFTOOLS_INDEX_2(BCFTOOLS_DEL_MLT_ALL.out.vcf)

    // Join biallelic VCF and TBI
    ch_biallelic_vcf_tbi = BCFTOOLS_DEL_MLT_ALL.out.vcf.join(BCFTOOLS_INDEX_2.out.tbi)

    // (Optional) Remove benchmarking samples (e.g. NA12878) from the reference panel
    if (!(params.remove_samples == null)){
        BCFTOOLS_REMOVE(ch_biallelic_vcf_tbi, [], [], [])
        BCFTOOLS_INDEX_3(BCFTOOLS_REMOVE.out.vcf)
        ch_biallelic_vcf_tbi = BCFTOOLS_REMOVE.out.vcf.join(BCFTOOLS_INDEX_3.out.tbi)
    }

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_biallelic_vcf_tbi, ch_fasta, [])

    // Output hap and legend files
    ch_hap_legend = BCFTOOLS_CONVERT.out.hap.join(BCFTOOLS_CONVERT.out.legend)

    emit:
    vcf_tbi        = ch_biallelic_vcf_tbi           // channel: [ [id, chr], vcf, tbi ]
    hap_legend     = ch_hap_legend                  // channel: [ [id, chr] '.hap', '.legend' ]
    versions       = ch_versions                    // channel: [ versions.yml ]
}