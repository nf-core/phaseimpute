workflow PREPARE_INPUT_STITCH {

    take:
    ch_input_impute  // channel:   [ [id, chr, region], bam, bai ]
    ch_posfile       // channel:   [ [panel, chr], sites, tsv ]
    ch_region        // channel:   [ [chr, region], region ]

    main:

    ch_versions      = Channel.empty()

    // Value channels
    def input_empty         = [[]]
    def rdata_empty         = [[]]
    k_val                   = params.k_val
    ngen                    = params.ngen

    // Get chromosomes of posfile
    ch_posfile = ch_posfile
        .map{metaPC, posfile -> [[chr: metaPC.chr], metaPC, posfile]}

    // Get chromosomes of fasta
    ch_chromosomes = ch_region
        .map{metaCR, region -> [[chr: metaCR.chr], metaCR.chr]}

    // Make final channel with parameters
    stitch_parameters = ch_posfile
        .map { it + input_empty + rdata_empty}
        .join(ch_chromosomes)
        .map { it + k_val + ngen}
        .map { metaC, metaPC, posfile, input, rdata, chr, k_val, ngen ->
            [metaPC, posfile, input, rdata, chr, k_val, ngen]
        }

    // Prepare sample files for STITCH
    // Group input by ID
    ch_bam_bai = ch_input_impute
        .map {metaI, bam, bai -> [metaI.subMap("id"), bam, bai]}
        .unique()

    // Make bamlist from bam input
    ch_bamlist = ch_bam_bai
        .map {it[1].toString().tokenize('/').last()}
        .collectFile(name: "bamlist.txt", newLine: true, sort: true)

    // Collect all files
    stitch_samples = ch_bam_bai
        .map {meta, bam, bai -> [[id: "all_samples"], bam, bai]}
        .groupTuple()
        .combine(ch_bamlist)
        .collect()

    emit:
    stitch_parameters = stitch_parameters // channel:   [ [chr], posfile, [], [], chr, k_val, ngen ]
    stitch_samples    = stitch_samples    // channel:   [ [id], bam, bai, bamlist ]
    versions          = ch_versions       // channel:   [ versions.yml ]
}
