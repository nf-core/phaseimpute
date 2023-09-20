include { BCFTOOLS_VIEW as VIEW_VCF_REGION       } from '../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_ANNOTATE                      } from '../../modules/nf-core/bcftools/annotate/main.nf'
include { BCFTOOLS_VIEW as VIEW_VCF_SNPS         } from '../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_VIEW as VIEW_VCF_SITES        } from '../../modules/nf-core/bcftools/view/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX1           } from '../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX2           } from '../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX3           } from '../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX4           } from '../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_INDEX5           } from '../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_NORM                          } from '../../modules/nf-core/bcftools/norm/main.nf'
include { BCFTOOLS_QUERY                         } from '../../modules/nf-core/bcftools/query/main.nf'
include { TABIX_BGZIP                            } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX                            } from '../../modules/nf-core/tabix/tabix/main'
include { VCF_PHASE_SHAPEIT5                     } from '../../subworkflows/nf-core/vcf_phase_shapeit5/main'


workflow GET_PANEL {
    take:
    ch_vcf          // channel: [ [id, ref], vcf ]
    ch_region       // channel: [ [ref, region], val(region) ]
    ch_fasta        // channel: [ fasta ]
    file_chr_rename // file

    main:

    ch_versions = Channel.empty()

    // Filter the region of interest of the panel file
    ch_input_region = ch_vcf
        .combine(ch_fasta)
        .combine(ch_region)
        .map{ metaI, vcf, index, fasta, metaR, region ->
            [metaI + metaR, vcf, index, region+",chr"+region]}

    VIEW_VCF_REGION(ch_input_region, [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_REGION.out.versions.first())

    VCF_INDEX1(VIEW_VCF_REGION.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX1.out.versions.first())

    // Rename the chromosome without prefix
    BCFTOOLS_ANNOTATE(VIEW_VCF_REGION.out.vcf
        .combine(VCF_INDEX1.out.csi, by:0)
        .combine(Channel.of([[],[], []])),
        file_chr_rename)
    
    VCF_INDEX2(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX2.out.versions.first())

    // Normalise the panel
    ch_norm = BCFTOOLS_ANNOTATE.out.vcf
        .combine(VCF_INDEX2.out.csi, by:0)
        .map{metaIRR, vcf, index -> [metaIRR.subMap(["ref","region"]), metaIRR, vcf, index]}
        .combine(ch_region, by:0)
        .map{metaRR, metaIRR, vcf, index, fasta, region ->
            [metaIRR, vcf, index, fasta]}

    BCFTOOLS_NORM(ch_norm)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    // Extract only the SNP
    VIEW_VCF_SNPS(BCFTOOLS_NORM.out.vcf
        .combine(Channel.of([[],[]])), [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SNPS.out.versions.first())

    VCF_INDEX3(VIEW_VCF_SNPS.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX3.out.versions.first())


    vcf_region = VIEW_VCF_SNPS.out.vcf
        .combine(VCF_INDEX3.out.csi, by:0)
    VIEW_VCF_SITES( vcf_region
        .combine(Channel.of([[]])),
        [], [], [])
    ch_versions = ch_versions.mix(VIEW_VCF_SITES.out.versions.first())

    VCF_INDEX4(VIEW_VCF_SITES.out.vcf)
    ch_versions = ch_versions.mix(VCF_INDEX4.out.versions.first())

    // Convert to TSV
    BCFTOOLS_QUERY(VIEW_VCF_SITES.out.vcf
        .combine(VCF_INDEX4.out.csi, by:0),
        [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    TABIX_BGZIP(BCFTOOLS_QUERY.out.txt)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    TABIX_TABIX(TABIX_BGZIP.out.output)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Phase panel
    vcf_region.view()
    VCF_PHASE_SHAPEIT5(vcf_region
            .map { meta, vcf, csi -> [meta, vcf, csi, [], meta.region] },
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[]]).collect())
    ch_versions = ch_versions.mix(VCF_PHASE_SHAPEIT5.out.versions.first())

    emit:
    panel_norm          = VIEW_VCF_SNPS.out.vcf
    panel_norm_index    = VCF_INDEX2.out.csi
    panel_sites         = VIEW_VCF_SITES.out.vcf
    panel_sites_index   = VCF_INDEX3.out.csi
    panel_tsv           = TABIX_BGZIP.out.output
    panel_tsv_index     = TABIX_TABIX.out.tbi
    panel_phased        = VCF_PHASE_SHAPEIT5.out.variants_phased
    panel_phased_index  = VCF_PHASE_SHAPEIT5.out.variants_index

    versions            = ch_versions      // channel: [ versions.yml ]
}
