include { BCFTOOLS_MPILEUP          } from '../../../modules/nf-core/bcftools/mpileup/main.nf'
include { BCFTOOLS_INDEX            } from '../../../modules/nf-core/bcftools/index/main.nf'


workflow COMPUTE_GL {

    take:
    ch_input   // channel: [ [id, ref], bam, bai ]
    ch_region  // channel: [ [ref, region], fasta, val(region)]
    ch_sites   // channel: [ [id, region], sites, index]
    ch_tsv     // channel: [ [id, region], tsv, index]

    main:

    ch_versions = Channel.empty()

    ch_panel = ch_sites
                .combine(ch_tsv, by:0)
    ch_mpileup = ch_input
        .map{ meta, bam, index -> [meta.subMap(["ref","region"]), meta, bam, index]}
        .combine(ch_region, by:0)
        .combine(ch_panel.map{metaIpRR,sites,tsv ->
                    [metaIpRR.subMap(["ref","region"]), metaIpRR, sites, tsv]},
                by:0)
        .map{metaRR, metaIRR, bam, bindex, fasta, region, metaIpRR, sites, tsv ->
                [metaIRR + ["panel": metaIpRR.panel], bam, fasta, sites, tsv]}

    BCFTOOLS_MPILEUP(ch_mpileup,[])
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    ch_output = BCFTOOLS_MPILEUP.out.vcf
        combine(BCFTOOLS_MPILEUP.out.tbi, by:0)    

    emit:
    vcf          = ch_output
    stats        = BCFTOOLS_MPILEUP.out.stats
    versions     = ch_versions      // channel: [ versions.yml ]
}