workflow CHANNEL_IMPUTE_CREATE_CSV {
    take:
    ch_impute_output     //  channel:   [ [id, tool], vcf, index ]
    outdir              //

    main:

    // Generate CSV from this step
    ch_impute_output.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/imputation/csv") { meta, vcf, index ->
        id              = meta.id
        vcf             = "${params.outdir}/imputation/${meta.tools}/concat/${vcf.fileName}" // Check if tool name is correct
        index           = "${params.outdir}/imputation/${meta.tools}/concat/${index.fileName}"
        tool            = meta.tools

        ["impute.csv", "id,vcf,index,tool\n${id},${vcf},${index},${tool}\n"]
    }


}
