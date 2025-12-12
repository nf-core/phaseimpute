// Helper functions for pipeline tests
class UTILS {
    public static def getPipelineResults(outdir, workflow){
        // stable_name: All files + folders in ${params.outdir}/ with a stable name
        def stable_name = getAllFilesFromDir(outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        // stable_path: All files in ${params.outdir}/ with stable content
        def stable_path = getAllFilesFromDir(outdir, ignoreFile: 'tests/.nftignore')
        // bam_files: All bam files
        def bam_files  = getAllFilesFromDir(outdir, include: ['**/*.bam'])
        // vcf_files: All vcf files
        def vcf_files  = getAllFilesFromDir(outdir, include: ['**/*.{vcf,bcf}.gz'])
        return [
            // Number of successful tasks
            workflow.trace.succeeded().size(),
            // pipeline versions.yml file for multiqc from which Nextflow version is removed because we tests pipelines on multiple Nextflow versions
            removeNextflowVersion("$outdir/pipeline_info/nf_core_phaseimpute_software_mqc_versions.yml"),
            // All stable path name, with a relative path
            stable_name,
            // All files with stable contents
            stable_path,
            // All bam files
            bam_files.collect { file -> [file.getName(), bam(file.toString()).readsMD5] },
            // All vcf files
            vcf_files.collect { file -> [
                file.getName(),
                path(file.toString()).vcf.variantsMD5
            ] }
        ]
    }
    public static def vcfDetails(filePath) {
        def summary = path(filePath).vcf.summary.replaceAll(", phasedAutodetect=(false|true)", "")
        def samples = path(filePath).vcf.header.getGenotypeSamples().sort()
        return [summary: summary, samples: samples]
    }
}
