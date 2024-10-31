include { getFileExtension         } from './subworkflows/local/utils_nfcore_phaseimpute_pipeline'

workflow {
    ext_array = getFileExtension([ file(params.pipelines_testdata_base_path + "hum_data/panel/chr21/1000GP.chr21.s.norel.vcf.gz", checkIfExists: true)])
    ext_file = getFileExtension(file(params.pipelines_testdata_base_path + "hum_data/panel/chr21/1000GP.chr21.s.norel.vcf.gz", checkIfExists: true))
    ext_empty = getFileExtension([ ])
    ext_string = getFileExtension(params.pipelines_testdata_base_path + "hum_data/panel/chr21/1000GP.chr21.s.norel.vcf.gz")
    ext_array2 = getFileExtension([params.pipelines_testdata_base_path + "hum_data/panel/chr21/1000GP.chr21.s.norel.vcf.gz", params.pipelines_testdata_base_path + "hum_data/panel/chr21/1000GP.chr21.s.norel.vcf.gz"])

    println ext_array
    println ext_file
    println ext_empty
    println ext_string
    println ext_array2
}
