include { BCFTOOLS_VIEW as VIEW_VCF_SNPS         } from '../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_VIEW as VIEW_VCF_SITES        } from '../../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX1           } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX3           } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX4           } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX5           } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_NORM                          } from '../../../modules/nf-core/bcftools/norm/main.nf'
include { BCFTOOLS_QUERY                         } from '../../../modules/nf-core/bcftools/query/main.nf'
include { TABIX_BGZIP                            } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX                            } from '../../../modules/nf-core/tabix/tabix/main'
include { VCF_PHASE_SHAPEIT5                     } from '../../../subworkflows/nf-core/vcf_phase_shapeit5/main'


workflow GET_PANEL {
    take:
    ch_vcf          // channel: [ [id], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_NORM(ch_vcf, ch_fasta.map{ genome, fasta, fai -> [genome, fasta] })
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    // Extract only the SNP
    VIEW_VCF_SNPS(BCFTOOLS_NORM.out.vcf // [ meta, vcf ]
        .combine(Channel.of([[]])), [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SNPS.out.versions.first())

    VCF_INDEX3(VIEW_VCF_SNPS.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX3.out.versions.first())

    ch_panel_norm = VIEW_VCF_SNPS.out.vcf
        .combine(VCF_INDEX3.out.csi, by:0)

    // Extract sites positions
    vcf_region = VIEW_VCF_SNPS.out.vcf
        .combine(VCF_INDEX3.out.csi, by:0)
    VIEW_VCF_SITES( ch_panel_norm,
        [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SITES.out.versions.first())

    VCF_INDEX4(VIEW_VCF_SITES.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX4.out.versions.first())

    ch_panel_sites = VIEW_VCF_SITES.out.vcf
        .combine(VCF_INDEX4.out.csi, by:0)

    // Convert to TSV
    BCFTOOLS_QUERY(ch_panel_sites,
        [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    TABIX_BGZIP(BCFTOOLS_QUERY.out.output)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    TABIX_TABIX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_panel_tsv = TABIX_BGZIP.out.output
        .combine(TABIX_TABIX.out.tbi, by: 0)

    // Phase panel
    if (params.phased == false) {
        VCF_PHASE_SHAPEIT5(vcf_region
            .map { meta, vcf, csi -> [meta, vcf, csi, [], meta.region] },
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[]]).collect())
        ch_versions = ch_versions.mix(VCF_PHASE_SHAPEIT5.out.versions.first())
        ch_panel_phased = VCF_PHASE_SHAPEIT5.out.variants_phased
            .combine(VCF_PHASE_SHAPEIT5.out.variants_index, by: 0)
    } else {
        ch_panel_phased = VIEW_VCF_SNPS.out.vcf
            .combine(VCF_INDEX3.out.csi, by: 0)
    }

    ch_panel = ch_panel_norm
        .combine(ch_panel_sites, by: 0)
        .combine(ch_panel_tsv, by: 0)
        .combine(ch_panel_phased, by: 0)
        .map{ metaI, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
            -> [[panel:metaI.id], norm, n_index, sites, s_index, tsv, t_index, phased, p_index]
        }

    emit:
    panel          = ch_panel         // channel: [ [panel], norm, n_index, sites, s_index, tsv, t_index, phased, p_index]
    versions       = ch_versions      // channel: [ versions.yml ]
}
