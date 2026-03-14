// Helper functions for pipeline tests
class UTILS {
    public static def getPipelineResults(outdir, workflow){
        // stable_name: All files + folders in ${params.outdir}/ with a stable name
        def stable_name = getAllFilesFromDir(outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        // stable_path: All files in ${params.outdir}/ with stable content
        def stable_path = getAllFilesFromDir(outdir, ignoreFile: 'tests/.nftignore')
        // bam_files: All bam files
        def bam_files   = getAllFilesFromDir(outdir, include: ['**/*.bam'])
        // vcf_files: All vcf files
        def vcf_files   = getAllFilesFromDir(outdir, include: ['**/*.{vcf,bcf}.gz'])
        // csv_files: All csv files
        def csv_files   = getAllFilesFromDir(outdir, include: ['**/*.csv'])
        return [
            // Number of successful tasks
            "workflow size": workflow.trace.succeeded().size(),
            // pipeline versions.yml file for multiqc from which Nextflow version is removed because we tests pipelines on multiple Nextflow versions
            "versions" : removeNextflowVersion("$outdir/pipeline_info/nf_core_phaseimpute_software_mqc_versions.yml"),
            // All stable path name, with a relative path
            "stable name": stable_name,
            // All files with stable contents
            "stable path": stable_path,
            // All bam files
            "BAM files": bam_files.collect { file -> [
                file.getName(),
                bam(file.toString()).readsMD5
            ] },
            // All vcf files
            "VCF files": vcf_files.collect { file -> [
                file.getName(),
                path(file.toString()).vcf.variantsMD5
            ] },
            // All csv files
            "CSV files": csv_files.collect { file ->
                def normalizedContent = path(file.toString()).text.replaceAll(/[^,\n]*\/([^,\n]+)/, '$1')
                [
                    fileName: file.getName(),
                    rows: normalizedContent.readLines()
                ]
            }
        ]
    }
    public static def vcfDetails(filePath) {
        def summary = path(filePath).vcf.summary.replaceAll(", phasedAutodetect=(false|true)", "")
        def samples = path(filePath).vcf.header.getGenotypeSamples().sort()
        return [summary: summary, samples: samples]
    }
}
