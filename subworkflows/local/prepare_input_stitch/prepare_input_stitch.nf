include { BCFTOOLS_QUERY                    } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_NORM                     } from '../../../modules/nf-core/bcftools/norm/main'
include { GAWK                              } from '../../../modules/nf-core/gawk'


workflow PREPARE_INPUT_STITCH {

    take:
    ch_panel_sites
    ch_fasta
    ch_input_impute

    main:

    ch_versions      = Channel.empty()

    // Prepare posfile and parameters for STITCH
    // Convert position file to tab-separated file
    BCFTOOLS_QUERY(ch_panel_sites, [], [], [])
    ch_posfile = BCFTOOLS_QUERY.out.output

    // Remove multiallelic positions from tsv
    GAWK(ch_posfile, [])

    // Get chromosomes of posfile
    ch_posfile = GAWK.out.output.map{meta, posfile -> return[['chr': meta.chr], posfile]}

    // Get chromosomes of fasta
    ch_chromosomes = ch_fasta.map{it -> it[2]}
                    .splitCsv(header: ["chr", "size", "offset", "lidebase", "linewidth", "qualoffset"], sep: "\t")
                    .map{it -> return [[chr: it.chr], it.chr]}

    // Combine channels
    def input_empty         = [[]]
    def rdata_empty         = [[]]
    k_val                   = params.k_val
    ngen                    = params.ngen

    // Make final channel with parameters
    stitch_parameters = ch_posfile.map { it + input_empty + rdata_empty}
                    .join(ch_chromosomes)
                    .map { it + k_val + ngen}

    // Prepare sample files for STITCH
    // Group input by ID
    ch_bam_bai = ch_input_impute.map {meta, bam, bai -> [[meta.id], bam, bai]}.unique()

    // Make bamlist from bam input
    ch_bamlist = ch_bam_bai
                    .map {it[1].tokenize('/').last()}
                    .collectFile(name: "bamlist.txt", newLine: true, sort: true)

    // Collect all files
    stitch_samples = ch_bam_bai.map {meta, bam, bai -> [["id": "all_samples"], bam, bai]}
                    .groupTuple()
                    .combine(ch_bamlist)
                    .collect()

    emit:
    stitch_parameters
    stitch_samples
    versions                   = ch_versions                           // channel:   [ versions.yml ]

}
