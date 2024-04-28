include { BCFTOOLS_VIEW as VIEW_VCF_SNPS         } from '../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_VIEW as VIEW_VCF_SITES        } from '../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX                         } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2     } from '../../../modules/nf-core/bcftools/index/main.nf'
include { TABIX_BGZIP                            } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX                            } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_QUERY                         } from '../../../modules/nf-core/bcftools/query/main.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_STITCH} from '../../../modules/nf-core/bcftools/query/main.nf'
include { GAWK as GAWK_STITCH                    } from '../../../modules/nf-core/gawk'



workflow VCF_SITES_EXTRACT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]

    main:

    ch_versions = Channel.empty()

    // Extract only SNPs from VCF
    VIEW_VCF_SNPS(ch_vcf, [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SNPS.out.versions.first())

    // Index SNPs
    BCFTOOLS_INDEX(VIEW_VCF_SNPS.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    // Join VCF and Index
    ch_panel_norm = VIEW_VCF_SNPS.out.vcf.combine(BCFTOOLS_INDEX.out.csi, by:0)

    // Extract sites positions
    VIEW_VCF_SITES( ch_panel_norm,[], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SITES.out.versions.first())

    // Index extracted sites
    BCFTOOLS_INDEX_2(VIEW_VCF_SITES.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())

    // Join extracted sites and index
    ch_panel_sites = VIEW_VCF_SITES.out.vcf.combine(BCFTOOLS_INDEX_2.out.csi, by:0)

    // Create empty channel

    ch_panel_tsv = []

    // Create TSVs for different tools

        // Convert to TSV with structure for Glimpse
        BCFTOOLS_QUERY(ch_panel_sites, [], [], [])
        ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

        // Compress TSV
        TABIX_BGZIP(BCFTOOLS_QUERY.out.output)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

        // Index compressed TSV
        TABIX_TABIX(TABIX_BGZIP.out.output)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        // Join compressed TSV and index
        ch_panel_tsv = TABIX_BGZIP.out.output.combine(TABIX_TABIX.out.tbi, by: 0)

        // TSV for STITCH
        // Convert position file to tab-separated file
        BCFTOOLS_QUERY_STITCH(ch_panel_sites, [], [], [])
        ch_posfile = BCFTOOLS_QUERY_STITCH.out.output

        // Remove multiallelic positions from tsv
        GAWK_STITCH(ch_posfile, [])

    emit:
    panel_tsv      = ch_panel_tsv
    vcf_tbi        = ch_panel_norm
    panel_sites    = ch_panel_sites
    posfile        = GAWK_STITCH.out.output
    versions       = ch_versions      // channel: [ versions.yml ]
}
