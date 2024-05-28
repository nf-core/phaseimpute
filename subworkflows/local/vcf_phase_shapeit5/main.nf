include { BEDTOOLS_MAKEWINDOWS                   } from '../../../modules/nf-core/bedtools/makewindows/main.nf'
include { SHAPEIT5_PHASECOMMON                   } from '../../../modules/nf-core/shapeit5/phasecommon/main'
include { SHAPEIT5_LIGATE                        } from '../../../modules/nf-core/shapeit5/ligate/main'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index/main.nf'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf        // channel (mandatory): [ [id, chr], vcf, csi, pedigree ]
    ch_region     // channel (optional) : [ [chr, region], region ]
    ch_ref        // channel (optional) : [ [id, chr], ref, csi ]
    ch_scaffold   // channel (optional) : [ [id, chr], scaffold, csi ]
    ch_map        // channel (optional) : [ [id, chr], map]

    main:

    ch_versions = Channel.empty()

    // It is needed to generate a file containing the region to phase in a Chr \tab Start \tab End format

    // Create the File in bed format and use the meta id for the file name
    ch_region_file = ch_region
        .collectFile(newLine: true) { metaCR, region -> ["${metaCR.chr}.bed", region.replace(":","\t").replace("-","\t")]}
        .map { file -> [[id: file.getBaseName(), chr:file.getBaseName()], file] }

    BEDTOOLS_MAKEWINDOWS(ch_region_file)
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.first())

    ch_chunk_output = BEDTOOLS_MAKEWINDOWS.out.bed
        .splitCsv(header: ['Chr', 'Start', 'End'], sep: "\t", skip: 0)
        .map { meta, it -> [meta.subMap("chr"), it["Chr"]+":"+it["Start"]+"-"+it["End"]]}

    ch_chunks_number = BEDTOOLS_MAKEWINDOWS.out.bed
        .map { meta, bed -> [meta.subMap("chr"), bed.countLines().intValue()]}

    ch_phase_input = ch_vcf
        .map { metaIC, vcf, index, pedigree ->
            [metaIC.subMap("chr"), metaIC, vcf, index, pedigree] }
        .combine(ch_chunk_output, by:0)
        .map { metaC, meta, vcf, index, pedigree, chunk ->
            [meta + [chunk: chunk], vcf, index, pedigree, chunk]
        }

    SHAPEIT5_PHASECOMMON (
        ch_phase_input, ch_ref,
        ch_scaffold, ch_map
    )
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_BCFTOOLS_INDEX_1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX_1.out.versions.first())

    ch_ligate_input = SHAPEIT5_PHASECOMMON.out.phased_variant
        .join(VCF_BCFTOOLS_INDEX_1.out.csi, failOnMismatch:true, failOnDuplicate:true)
        .map{ meta, vcf, csi -> [meta.subMap("id", "chr"), [vcf, meta.chunk], csi]}
        .groupTuple()
        .view()
        .map{ meta, vcf, csi ->
                [ meta,
                vcf
                    .sort { a, b ->
                        def aStart = a.last().split("-")[-1].toInteger()
                        def bStart = b.last().split("-")[-1].toInteger()
                        aStart <=> bStart
                    }
                    .collect{it.first()},
                csi]}
        .view()

    SHAPEIT5_LIGATE(ch_ligate_input)
    ch_versions = ch_versions.mix(SHAPEIT5_LIGATE.out.versions.first())

    VCF_BCFTOOLS_INDEX_2(SHAPEIT5_LIGATE.out.merged_variants)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX_2.out.versions.first())

    ch_vcf_tbi_join = SHAPEIT5_LIGATE.out.merged_variants
        .join(VCF_BCFTOOLS_INDEX_2.out.csi)

    emit:
    vcf_tbi_join        = ch_vcf_tbi_join         // channel: [ [id, chr], vcf, csi ]
    versions            = ch_versions             // channel: [ versions.yml ]
}

