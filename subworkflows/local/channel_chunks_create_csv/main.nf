workflow CHANNEL_CHUNKS_CREATE_CSV {
    take:
    ch_chunks           //  channel:   [ [id, chr], tsv ]
    outdir              //

    main:

    // Generate CSV from this step
    ch_chunks.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/prep_panel/csv") { meta, file ->
        chr             =  meta.chr
        file            = "${params.outdir}/prep_panel/chunks/${file.fileName}"

        ["chunks.csv", "chr,file\n${chr},${file}\n"]
    }


}
