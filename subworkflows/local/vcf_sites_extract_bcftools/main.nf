include { BCFTOOLS_CONVERT              } from '../../../modules/nf-core/bcftools/convert'
include { BCFTOOLS_VIEW                 } from '../../../modules/nf-core/bcftools/view'

workflow VCF_SITES_EXTRACT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()
    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_vcf, ch_fasta, [])
    ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions.first())

    // Extract sites positions
    BCFTOOLS_VIEW(ch_vcf, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Join extracted sites and index
    ch_posfile = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi)
        .join(BCFTOOLS_CONVERT.out.hap)
        .join(BCFTOOLS_CONVERT.out.legend)

    emit:
    posfile       = ch_posfile          // channel: [ [id, chr], vcf, csi, hap, legend ]
    versions      = ch_versions         // channel: [ versions.yml ]
}
