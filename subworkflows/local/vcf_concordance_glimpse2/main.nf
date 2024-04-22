include { GLIMPSE2_CONCORDANCE        } from '../../../modules/nf-core/glimpse2/concordance'
include { GAWK as CONCATENATE         } from '../../../modules/nf-core/gawk'
include { ADD_COLUMNS                 } from '../../../modules/local/addcolumns'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip'

workflow VCF_CONCORDANCE_GLIMPSE2 {

    take:
        ch_vcf_emul   // VCF file with imputed genotypes [[id, chr, region, panel, simulate, tools], vcf, csi]
        ch_vcf_truth  // VCF file with truth genotypes   [[id, chr, region], vcf, csi]
        ch_vcf_freq   // VCF file with panel frequencies [[panel, chr], vcf, csi]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_concordance = ch_vcf_emul
        .map{
            metaICRPST, vcf, csi ->
            [metaICRPST.subMap(["id", "chr", "region"]), metaICRPST, vcf, csi]
        }
        .combine(ch_vcf_truth.map{
                metaICRP, vcf, csi ->
                [metaICRP.subMap(["id", "chr", "region"]), vcf, csi]
            },
            by:0
        )
        .map{metaICR, metaIPCRTS, emul, e_csi, truth, t_csi ->
            [metaICR.subMap(["chr"]), metaIPCRTS, emul, e_csi, truth, t_csi]
        }
        .combine(ch_vcf_freq.map{metaCRP, vcf, csi ->
                    [metaCRP.subMap(["chr"]), metaCRP, vcf, csi]},
                by:0)
        .map{metaC, metaIPCRTS, emul, e_csi, truth, t_csi, metaCRP, freq, f_csi ->
            [metaIPCRTS, emul, e_csi, truth, t_csi, freq, f_csi, [], metaIPCRTS.region]
        }

    GLIMPSE2_CONCORDANCE (
        ch_concordance,
        [[], [], params.bins, [], []],
        params.min_val_gl, params.min_val_dp
    )
    ch_versions = ch_versions.mix(GLIMPSE2_CONCORDANCE.out.versions.first())

    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_cal.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_grp.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_spl.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_grp.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_spl.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_per_site.map{meta, txt -> [txt]})

    GUNZIP(GLIMPSE2_CONCORDANCE.out.errors_grp)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())
    ADD_COLUMNS(GUNZIP.out.gunzip)
    ch_versions = ch_versions.mix(ADD_COLUMNS.out.versions.first())

    CONCATENATE(
        ADD_COLUMNS.out.txt
            .map{meta, txt -> [["id":"TestQuality"], txt]}
            .groupTuple(),
        Channel.of(
            '(NR == 1) || (FNR > 1)'
        ).collectFile(name:"program.txt")
    )
    ch_versions = ch_versions.mix(CONCATENATE.out.versions.first())

    emit:
    stats           = CONCATENATE.out.output      // [ meta, txt ]
    versions        = ch_versions                 // channel: [ versions.yml ]
    multiqc_files   = ch_multiqc_files
}
