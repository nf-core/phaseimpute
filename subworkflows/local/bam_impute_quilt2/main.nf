include { QUILT_QUILT2    } from '../../../modules/nf-core/quilt/quilt2'
include { GLIMPSE2_LIGATE } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX  } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_QUILT2 {
    take:
    ch_input // channel (mandatory):   [ [id], [bam], [bai], bampaths, bamnames ]
    ch_reference_panel // channel (mandatory):   [ [panel, chr], vcf, index ]
    ch_chunks // channel (optional) :   [ [panel, chr], chr, start, end ]
    ch_map // channel (optional) :   [ [panel, chr], map ]
    ch_fasta // channel (optional) :   [ [genome], fa, fai ]
    n_gen // integer: Number of generations since founding or mixing
    buffer // integer: Buffer of region to perform imputation over

    main:

    ch_parameters = ch_reference_panel
        .combine(ch_map, by: 0)
        .combine(ch_chunks, by: 0)

    ch_parameters.ifEmpty {
        error("ERROR: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input.")
    }

    ch_bam_params = ch_input
        .combine(ch_parameters)
        .map { metaI, bam, bai, bampath, bamname, metaPC, reference_vcf, reference_index, gmap, chr, start, end ->
            def regionout = "${chr}"
            if (start != [] && end != []) {
                regionout = "${chr}:${start}-${end}"
            }
            [
                metaPC + metaI + ["regionout": regionout],
                bam,
                bai,
                bampath,
                bamname,
                reference_vcf,
                reference_index,
                [],
                [],
                [],
                chr,
                start,
                end,
                n_gen,
                buffer,
                gmap,
            ]
        }

    QUILT_QUILT2(ch_bam_params, ch_fasta)

    ligate_input = QUILT_QUILT2.out.vcf
        .join(QUILT_QUILT2.out.tbi)
        .map { meta, vcf, index ->
            def keysToKeep = meta.keySet() - ['regionout']
            [meta.subMap(keysToKeep), vcf, index]
        }
        .groupTuple()

    GLIMPSE2_LIGATE(ligate_input)

    BCFTOOLS_INDEX(GLIMPSE2_LIGATE.out.merged_variants)

    ch_vcf_index = GLIMPSE2_LIGATE.out.merged_variants.join(
        BCFTOOLS_INDEX.out.tbi.mix(BCFTOOLS_INDEX.out.csi),
        failOnMismatch: true,
        failOnDuplicate: true,
    )

    emit:
    vcf_index = ch_vcf_index // channel:   [ [id, chr], vcf, tbi ]
}
