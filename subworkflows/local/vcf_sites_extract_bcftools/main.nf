include { BCFTOOLS_CONVERT              } from '../../../modules/nf-core/bcftools/convert'
include { BCFTOOLS_VIEW                 } from '../../../modules/nf-core/bcftools/view'
include { GAWK                          } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP                   } from '../../../modules/nf-core/tabix/bgzip'

workflow VCF_SITES_EXTRACT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_vcf, ch_fasta, [])

    // Extract sites positions
    BCFTOOLS_VIEW(ch_vcf, [], [], [])

    // Transform posfile to TSV with ','
    GAWK(BCFTOOLS_CONVERT.out.legend, [], false)

    // Compress TSV
    TABIX_BGZIP(GAWK.out.output)

    // Join extracted sites and index
    ch_posfile = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi)
        .join(BCFTOOLS_CONVERT.out.hap)
        .join(BCFTOOLS_CONVERT.out.legend)
        .join(TABIX_BGZIP.out.output)

    emit:
    posfile       = ch_posfile          // channel: [ [id, chr], vcf, csi, hap, legend, posfile ]
}
