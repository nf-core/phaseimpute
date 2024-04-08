include { BCFTOOLS_INDEX                    } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2} from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_3} from '../../../modules/nf-core/bcftools/index/main'
include { GLIMPSE_CHUNK                     } from '../../../modules/nf-core/glimpse/chunk/main'
include { BCFTOOLS_CONVERT                  } from '../../../modules/nf-core/bcftools/convert/main'
include { BCFTOOLS_NORM                     } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW                     } from '../../../modules/nf-core/bcftools/view/main'


workflow MAKE_CHUNKS {

    take:
    ch_reference                           // channel: [ val(meta),vcf ]

    main:

    ch_versions = Channel.empty()

    // Make chunks
    ch_vcf_csi_chr = ch_reference.map{meta, vcf, csi -> [meta, vcf, csi, meta.chr]}
    GLIMPSE_CHUNK(ch_vcf_csi_chr)

    // Rearrange chunks into channel
    ch_chunks = GLIMPSE_CHUNK.out.chunk_chr
                    .splitText()
                    .map { metamap, line ->
                        def fields = line.split("\t")
                        def startEnd = fields[2].split(':')[1].split('-')
                        [metamap, metamap.chr, startEnd[0], startEnd[1]]
                    }

    ch_fasta = Channel.value([file(params.fasta,checkIfExists: true)])
                            .map { file -> [[id: 'genome'], file] }

    // Join duplicated biallelic sites into multiallelic records
    BCFTOOLS_NORM(ch_reference, ch_fasta)

    // Index multiallelic VCF
    BCFTOOLS_INDEX_2(BCFTOOLS_NORM.out.vcf)

    // Join multiallelic VCF and TBI
    ch_multiallelic_vcf_tbi = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_INDEX_2.out.tbi)

    // Remove all multiallelic records:
    BCFTOOLS_VIEW(ch_multiallelic_vcf_tbi, [], [], [])

    // Index biallelic VCF
    BCFTOOLS_INDEX_3(BCFTOOLS_VIEW.out.vcf)

    // Join biallelic VCF and TBI
    ch_biallelic_vcf_tbi = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_INDEX_3.out.tbi)

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_biallelic_vcf_tbi, ch_fasta, [])

    // Output hap and legend files
    ch_hap_legend = BCFTOOLS_CONVERT.out.hap.join(BCFTOOLS_CONVERT.out.legend)





    emit:
    ch_chunks                 = ch_chunks                              // channel: [ chr, val(meta), start, end, number ]
    ch_hap_legend             = ch_hap_legend                          // channel: [ chr, val(meta), hap, legend ]
    versions                  = ch_versions                           // channel:  [ versions.yml ]
}
