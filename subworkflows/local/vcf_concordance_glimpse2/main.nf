include { GLIMPSE2_CONCORDANCE        } from '../../../modules/nf-core/glimpse2/concordance'
include { GAWK                        } from '../../../modules/nf-core/gawk'
include { ADDCOLUMNS                 } from '../../../modules/local/addcolumns'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip'
include { GAWK as GAWK_ERROR_SPL      } from '../../../modules/nf-core/gawk'
include { GAWK as GAWK_RSQUARE_SPL    } from '../../../modules/nf-core/gawk'

workflow VCF_CONCORDANCE_GLIMPSE2 {

    take:
        ch_vcf_emul   // VCF file with imputed genotypes [ [id, panel, tool], vcf, csi]
        ch_vcf_truth  // VCF file with truth genotypes   [ [id, panel], vcf, csi]
        ch_vcf_freq   // VCF file with panel frequencies [ [panel_id, chr], vcf, csi]
        ch_region     // Regions to process              [ [chr, region], region]

    main:

    ch_multiqc_files = channel.empty()

    ch_concordance = ch_vcf_emul
        .map{metaIPTC, vcf, csi -> [metaIPTC.subMap("id"), metaIPTC, vcf, csi]}
        .combine(ch_vcf_truth
            .map{metaIPC, vcf, csi -> [ metaIPC.subMap("id"), vcf, csi ]}
            , by: 0
        )
        .combine(ch_vcf_freq)
        .combine(ch_region.map{ _meta, region ->
            [ region ]
        }.collect().toList())
        .map{ _metaI, metaIPTC, emul, e_csi, truth, t_csi, _metaP, freq, f_csi, regions ->
            [metaIPTC, emul, e_csi, truth, t_csi, freq, f_csi, [], regions]
        }

    GLIMPSE2_CONCORDANCE (
        ch_concordance,
        [[], [], params.bins, [], [], params.min_val_gl, params.min_val_dp]
    )

    GAWK_ERROR_SPL(
        GLIMPSE2_CONCORDANCE.out.errors_spl,
        [],
        false
    )

    GAWK_RSQUARE_SPL(
        GLIMPSE2_CONCORDANCE.out.rsquare_spl,
        [],
        false
    )

    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_cal.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_grp.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GAWK_ERROR_SPL.out.output.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_grp.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GAWK_RSQUARE_SPL.out.output.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_per_site.map{ _meta, txt -> [txt]})

    GUNZIP(GLIMPSE2_CONCORDANCE.out.errors_grp)

    ADDCOLUMNS(GUNZIP.out.gunzip)

    GAWK(
        ADDCOLUMNS.out.txt
            .map { meta, files ->
                // Normalize tools to always be a list
                def normalizedMeta = meta.clone()
                normalizedMeta.tools = (meta.tools instanceof List) ? meta.tools : [meta.tools]
                [normalizedMeta, files]
            }
            .toSortedList { a, b ->
                def keyA = "${a[0].id},${a[0].tools.sort().join(',')}"
                def keyB = "${b[0].id},${b[0].tools.sort().join(',')}"
                keyA <=> keyB
            }
            .map { sorted_list ->
                def all_files = sorted_list.collect { _meta, file -> file }
                [["id": "AllSamples"], all_files]
            },
        [],
        false
    )

    emit:
    stats           = GAWK.out.output             // [ [all], txt ]
    multiqc_files   = ch_multiqc_files
}
