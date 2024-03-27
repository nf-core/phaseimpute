include { GLIMPSE2_CONCORDANCE        } from '../../../modules/nf-core/glimpse2/concordance'
include { CONCATENATE                 } from '../../../modules/local/concatenate'
include { ADD_COLUMNS                 } from '../../../modules/local/addcolumns'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip'

workflow VCF_CONCORDANCE_GLIMPSE {

    take:
        ch_vcf_emul   // VCF file with imputed genotypes [[id, chr, region, panel, simulate, tools], vcf, csi]
        ch_vcf_truth  // VCF file with truth genotypes   [[id, chr, region], vcf, csi]
        ch_vcf_freq   // VCF file with panel frequencies [[panel, chr], vcf, csi]

    main:

    ch_versions = Channel.empty()

    ch_concordance = ch_vcf_emul
        .map{
            metaICRPST, vcf, csi ->
            [metaICRPST.subMap(["id", "chr", "region", "panel"]), metaICRPST, vcf, csi]
        }
        .combine(ch_vcf_truth, by:0)
        .map{metaICRP, metaIPCRTS, emul, e_csi, truth, t_csi ->
            [metaICRP.subMap(["chr"]), metaIPCRTS, emul, e_csi, truth, t_csi]
        }
        .combine(ch_vcf_freq.map{metaCRP, vcf, csi ->
                    [metaCRP.subMap(["chr"]), metaCRP, vcf, csi]},
                by:0)
        .map{metaC, metaIPCRTS, emul, e_csi, truth, t_csi, metaCRP, freq, f_csi ->
            [metaIPCRTS, emul, e_csi, truth, t_csi, freq, f_csi, [], metaIPCRTS.region]
        }

    GLIMPSE2_CONCORDANCE (
        ch_concordance,
        [[], [], "0 0.01 0.05 0.1 0.2 0.5", [], []],
        0.9, 5
    )
    GUNZIP(GLIMPSE2_CONCORDANCE.out.errors_grp)
    ADD_COLUMNS(GUNZIP.out.gunzip)

    CONCATENATE(ADD_COLUMNS.out.txt
                    .map{meta, txt -> [["id":"TestQuality"], txt]}
                    .groupTuple()
    )


    emit:
    stats     = CONCATENATE.out.txt         // [ meta, txt ]
    versions  = ch_versions                 // channel: [ versions.yml ]
}
