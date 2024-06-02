workflow CHANNEL_SIMULATE_CREATE_CSV {
    take:
    ch_input_impute     //  channel: [ [id, genome, chr, region, depth], bam, bai ]
    outdir              //

    main:
    // Generate CSV from this step
    ch_input_impute.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/simulation/csv") { meta, file, index ->
        sample          =  meta.id
        file            = "${params.outdir}/simulation/${file.fileName}"
        index           = "${params.outdir}/simulation/${index.fileName}"

        ["simulate.csv", "sample,file,index\n${sample},${file},${index}\n"]
    }


}
