include { BCFTOOLS_VIEW                 } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_INDEX                } from '../../../modules/nf-core/bcftools/index'
include { TABIX_BGZIP                   } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_TABIX                   } from '../../../modules/nf-core/tabix/tabix'
include { BCFTOOLS_QUERY                } from '../../../modules/nf-core/bcftools/query'
include { GAWK                          } from '../../../modules/nf-core/gawk'
include { BCFTOOLS_CONVERT              } from '../../../modules/nf-core/bcftools/convert'


workflow VCF_SITES_EXTRACT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()
    ch_fasta = ch_fasta.map { meta, fasta, fai -> [meta, fasta] }

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_vcf, ch_fasta, [])
    ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions)

    // Output hap and legend files
    ch_hap_legend = BCFTOOLS_CONVERT.out.hap.join(BCFTOOLS_CONVERT.out.legend)

    // Extract sites positions
    BCFTOOLS_VIEW(ch_vcf, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Index extracted sites
    BCFTOOLS_INDEX(BCFTOOLS_VIEW.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    // Join extracted sites and index
    ch_panel_sites = BCFTOOLS_VIEW.out.vcf.combine(BCFTOOLS_INDEX.out.csi, by:0)

    // Convert to TSV with structure for Glimpse
    BCFTOOLS_QUERY(ch_panel_sites, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    // Compress TSV
    TABIX_BGZIP(BCFTOOLS_QUERY.out.output)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    // Generate default posfile (sites vcf, sites index and sites txt)
    ch_posfile = ch_panel_sites
            .join(TABIX_BGZIP.out.output)

    // Convert TSV to Stitch format ","" to "\t"
    GAWK(BCFTOOLS_QUERY.out.output, [])
    ch_versions = ch_versions.mix(GAWK.out.versions)

    // Generate glimpse posfile
    ch_glimpse_posfile = ch_posfile.map{ metaPC, sites, s_index, tsv -> [metaPC, sites, tsv]}

    emit:
    hap_legend             = ch_hap_legend       // channel: [ [id, chr], '.hap', '.legend' ]
    panel_tsv_stitch       = GAWK.out.output     // channel: [ [id, chr], txt ]
    panel_sites            = ch_panel_sites      // channel: [ [id, chr], vcf, csi ]
    posfile                = ch_posfile          // channel: [ [id, chr], vcf, csi, tsv.gz ]
    glimpse_posfile        = ch_glimpse_posfile  // channel: [ [id, chr], vcf, tsv.gz ]
    versions               = ch_versions         // channel: [ versions.yml ]
}
