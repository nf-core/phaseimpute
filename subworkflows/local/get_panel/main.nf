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
    ch_vcf          // channel: [ [id, ref], vcf, index ]
    ch_region       // channel: [ [ref, region], val(region) ]
    ch_fasta        // channel: [ fasta ]
    file_chr_rename // file rename

    main:

    ch_versions = Channel.empty()

    // Normalise the panel
    ch_norm = ch_vcf
        .combine(ch_fasta)

    BCFTOOLS_NORM(ch_norm)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    // Extract only the SNP
    VIEW_VCF_SNPS(BCFTOOLS_NORM.out.vcf
        .combine(Channel.of([[],[]])), [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SNPS.out.versions.first())

    VCF_INDEX3(VIEW_VCF_SNPS.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX3.out.versions.first())

    ch_panel_norm = VIEW_VCF_SNPS.out.vcf
        .combine(VCF_INDEX3, by:0)

    // Extract sites positions
    vcf_region = VIEW_VCF_SNPS.out.vcf
        .combine(VCF_INDEX3.out.csi, by:0)
    VIEW_VCF_SITES( vcf_region
        .combine(Channel.of([[]])),
        [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SITES.out.versions.first())

    VCF_INDEX4(VIEW_VCF_SITES.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX4.out.versions.first())

    ch_panel_sites = VIEW_VCF_SITES.out.vcf
        .combine(VCF_INDEX4, by:0)

    // Convert to TSV
    BCFTOOLS_QUERY(VIEW_VCF_SITES.out.vcf
        .combine(VCF_INDEX4.out.csi, by:0),
        [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    TABIX_BGZIP(BCFTOOLS_QUERY.out.txt)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    TABIX_TABIX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_panel_tsv = TABIX_BGZIP.out.output
        .combine(TABIX_TABIX.out.tbi, by: 0)

    // Phase panel
    if (params.phase_panel == true) {
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

    emit:
    panel_norm          = ch_panel_norm    // channel: [ meta, vcf, index ]
    panel_sites         = ch_panel_sites   // channel: [ meta, bcf, index ]
    panel_tsv           = ch_panel_tsv     // channel: [ meta, tsv, index ]
    panel_phased        = ch_panel_phased  // channel: [ meta, vcf, index ]

    versions            = ch_versions      // channel: [ versions.yml ]
}
