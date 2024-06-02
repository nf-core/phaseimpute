workflow CHANNEL_PANEL_CREATE_CSV {
    take:
    ch_panel           //  channel:   [ [panel, chr], vcf, index ]
    ch_hap_legend      //  channel:   [ [id, chr], '.hap', '.legend' ]
    outdir             //

    main:

    // Generate CSV from this step
    ch_panel.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/prep_panel/csv") { meta, vcf, index ->
        panel           = meta.panel
        chr             = meta.chr
        vcf             = "${params.outdir}/prep_panel/panel/${vcf.fileName}" // This is currently not being exported
        index           = "${params.outdir}/prep_panel/panel/${index.fileName}"


        ["panel.csv", "panel,chr,vcf,index\n${panel},${chr},${vcf},${index}\n"]
    }


}
