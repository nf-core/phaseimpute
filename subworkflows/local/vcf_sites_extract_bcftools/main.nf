include { BCFTOOLS_VIEW as VIEW_VCF_SITES        } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2     } from '../../../modules/nf-core/bcftools/index'
include { TABIX_BGZIP                            } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_TABIX                            } from '../../../modules/nf-core/tabix/tabix'
include { BCFTOOLS_QUERY                         } from '../../../modules/nf-core/bcftools/query'

workflow VCF_SITES_EXTRACT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]

    main:

    ch_versions = Channel.empty()

    // Extract sites positions
    VIEW_VCF_SITES( ch_vcf,[], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SITES.out.versions.first())

    // Index extracted sites
    BCFTOOLS_INDEX_2(VIEW_VCF_SITES.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())

    // Join extracted sites and index
    ch_panel_sites = VIEW_VCF_SITES.out.vcf.combine(BCFTOOLS_INDEX_2.out.csi, by:0)

    // Convert to TSV with structure for Glimpse
    BCFTOOLS_QUERY_TSV(ch_panel_sites, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY_TSV.out.versions.first())

    // Compress TSV
    TABIX_BGZIP(BCFTOOLS_QUERY_TSV.out.output)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    // Index compressed TSV
    TABIX_TABIX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Join compressed TSV and index
    ch_panel_tsv = TABIX_BGZIP.out.output.combine(TABIX_TABIX.out.tbi, by: 0)

    emit:
    panel_tsv      = ch_panel_tsv
    vcf_tbi        = ch_vcf
    panel_sites    = ch_panel_sites
    versions       = ch_versions      // channel: [ versions.yml ]
}
