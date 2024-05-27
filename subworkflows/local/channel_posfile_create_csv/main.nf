workflow CHANNEL_POSFILE_CREATE_CSV {
    take:
    ch_posfile           //  channel:   [ [id, chr], tsv ]
    outdir              //

    main:
    // Generate CSV from this step
    ch_posfile.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/prep_panel/csv") { meta, file ->
        chr             =  meta.chr
        file            = "${params.outdir}/prep_panel/posfile/${file.fileName}"

        ["posfile.csv", "chr,file\n${chr},${file}\n"]
    }


}
