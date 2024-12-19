//
// Subworkflow with functionality specific to the nf-core/phaseimpute pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { GET_REGION                } from '../get_region'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions      = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )


    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )


    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create fasta channel
    //

    genome = params.genome ? params.genome : file(params.fasta, checkIfExists:true).getBaseName()
    if (params.genome) {
        genome = params.genome
        ch_fasta  = Channel.of([[genome:genome], getGenomeAttribute('fasta')])
        fai       = getGenomeAttribute('fai')
        if (fai == null) {
            SAMTOOLS_FAIDX(ch_fasta, Channel.of([[], []]))
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
            fai         = SAMTOOLS_FAIDX.out.fai.map{ it[1] }
        } else {
            fai = Channel.of(file(fai, checkIfExists:true))
        }
    } else if (params.fasta) {
        genome = file(params.fasta, checkIfExists:true).getBaseName()
        ch_fasta  = Channel.of([[genome:genome], file(params.fasta, checkIfExists:true)])
        if (params.fasta_fai) {
            fai = Channel.of(file(params.fasta_fai, checkIfExists:true))
        } else {
            SAMTOOLS_FAIDX(ch_fasta, Channel.of([[], []]))
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
            fai         = SAMTOOLS_FAIDX.out.fai.map{ it[1] }
        }
    }
    ch_ref_gen = ch_fasta.combine(fai).collect()

    //
    // Create channel from input file provided through params.input
    //
    if (params.input) {
        ch_input = Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map { meta, file, index -> [meta + [batch: 0], file, index] } // Set batch to 0 by default
    } else {
        ch_input = Channel.of([[], [], []])
    }

    // Check that the batch size and extension is compatible with the tools
    validateInputBatchTools(
        ch_input,
        params.batch_size,
        getFilesSameExt(ch_input),
        params.tools ? params.tools.split(',') : []
    )

    //
    // Create channel from input file provided through params.input_truth
    //
    if (params.input_truth) {
        if (params.input_truth.endsWith("csv")) {
            ch_input_truth = Channel
                .fromList(samplesheetToList(params.input_truth, "${projectDir}/assets/schema_input.json"))
                .map {
                    meta, file, index ->
                        [ meta, file, index ]
                }
            // Check if all extension are identical
            getFilesSameExt(ch_input_truth)
        } else {
            // #TODO Wait for `oneOf()` to be supported in the nextflow_schema.json
            error "Panel file provided is of another format than CSV (not yet supported). Please separate your panel by chromosome and use the samplesheet format."
        }
    } else {
        ch_input_truth = Channel.of([[], [], []])
    }

    //
    // Create channel for panel
    //
    if (params.panel) {
        if (params.panel.endsWith("csv")) {
            println "Panel file provided as input is a samplesheet"
            ch_panel = Channel.fromList(samplesheetToList(
                params.panel, "${projectDir}/assets/schema_input_panel.json"
            ))
        } else {
            // #TODO Wait for `oneOf()` to be supported in the nextflow_schema.json
            error "Panel file provided is of another format than CSV (not yet supported). Please separate your panel by chromosome and use the samplesheet format."
        }
    } else {
        // #TODO check if panel is required
        ch_panel        = Channel.of([[],[],[]])
    }

    //
    // Create channel from region input
    //
    if (params.input_region == null){
        // #TODO Add support for string input
        GET_REGION ("all", ch_ref_gen)
        ch_versions = ch_versions.mix(GET_REGION.out.versions)
        ch_regions  = GET_REGION.out.regions
    }  else  if (params.input_region.endsWith(".csv")) {
        println "Region file provided as input is a csv file"
        ch_regions = Channel.fromList(samplesheetToList(
            params.input_region, "${projectDir}/assets/schema_input_region.json"
        ))
        .map{ chr, start, end ->
            assert end >= start : "End position must be greater than or equal to start position"
            [["chr": chr], chr + ":" + start + "-" + end]
        }
        .map{ metaC, region -> [metaC + ["region": region], region]}
    } else {
        error "Region file provided is of another format than CSV (not yet supported). Please separate your reference genome by chromosome and use the samplesheet format."
    }

    //
    // Create map channel
    //
    if (params.map) {
        if (params.map.endsWith(".csv")) {
            println "Map file provided as input is a samplesheet"
            ch_map = Channel.fromList(samplesheetToList(params.map, "${projectDir}/assets/schema_map.json"))
        } else {
            error "Map file provided is of another format than CSV (not yet supported). Please separate your reference genome by chromosome and use the samplesheet format."
        }
    } else {
        ch_map = ch_regions
            .map{ metaCR, _regions -> [metaCR.subMap("chr"), []] }
    }

    //
    // Create depth channel
    //
    if (params.depth) {
        ch_depth = Channel.of([[depth: params.depth], params.depth])
    } else {
        ch_depth = Channel.of([[],[]])
    }

    //
    // Create genotype array channel
    //
    if (params.genotype) {
        ch_genotype = Channel.of([[gparray: params.genotype], params.genotype])
    } else {
        ch_genotype = Channel.of([[],[]])
    }

    //
    // Create posfile channel
    //
    if (params.posfile) {
        ch_posfile = Channel // ["panel", "chr", "vcf", "index", "hap", "legend"]
            .fromList(samplesheetToList(params.posfile, "${projectDir}/assets/schema_posfile.json"))
    } else {
        ch_posfile = Channel.of([[],[],[],[],[]])
    }

    if (!params.steps.split(',').contains("panelprep") & !params.steps.split(',').contains("all")) {
        validatePosfileTools(
            ch_posfile,
            params.tools ? params.tools.split(','): [],
            params.steps.split(',')
        )
    }

    //
    // Create chunks channel
    //
    if (params.chunks) {
        ch_chunks = Channel
            .fromList(samplesheetToList(params.chunks, "${projectDir}/assets/schema_chunks.json"))
    } else {
        ch_chunks = Channel.of([[],[]])
    }

    //
    // Check contigs name in different meta map
    //
    // Collect all chromosomes names in all different inputs
    chr_ref = ch_ref_gen.map { _meta, _fasta, fai_file -> [fai_file.readLines()*.split('\t').collect{it[0]}] }
    chr_regions = extractChr(ch_regions)

    // Check that the chromosomes names that will be used are all present in different inputs
    chr_ref_mis     = checkMetaChr(chr_regions, chr_ref, "reference genome")
    chr_chunks_mis  = checkMetaChr(chr_regions, extractChr(ch_chunks), "chromosome chunks")
    chr_map_mis     = checkMetaChr(chr_regions, extractChr(ch_map), "genetic map")
    chr_panel_mis   = checkMetaChr(chr_regions, extractChr(ch_panel), "reference panel")
    chr_posfile_mis = checkMetaChr(chr_regions, extractChr(ch_posfile), "position")

    // Compute the intersection of all chromosomes names
    chr_all_mis = chr_ref_mis.concat(chr_chunks_mis, chr_map_mis, chr_panel_mis, chr_posfile_mis)
        .unique()
        .toList()
        .subscribe{ chr ->
            if (chr.size() > 0) {
                def chr_names = chr.size() > params.max_chr_names ? chr[0..params.max_chr_names - 1] + ['...'] : chr
                log.warn "The following contigs are absent from at least one file : ${chr_names} and therefore won't be used" } }

    ch_regions = ch_regions
        .combine(chr_all_mis.toList())
        .filter { meta, _regions, chr_mis ->
            !(meta.chr in chr_mis)
        }
        .map { meta, regions, _chr_mis -> [meta, regions] }
        .ifEmpty { error "No regions left to process" }

    ch_regions
        .map { it[1] }
        .collect()
        .subscribe { log.info "The following contigs will be processed: ${it}" }

    // Remove other contigs from panel and posfile files
    ch_panel = ch_panel
        .combine(ch_regions.collect{ it[0]["chr"]}.toList())
        .filter { meta, _vcf, _index, chrs ->
            meta.chr in chrs
        }
        .map {meta, vcf, index, _chrs ->
            [meta, vcf, index]
        }

    ch_posfile = ch_posfile
        .combine(ch_regions.collect{ it[0]["chr"]}.toList())
        .filter { meta, _vcf, _index, _hap, _legend, chrs ->
            meta.chr in chrs
        }
        .map {meta, vcf, index, hap, legend, _chrs ->
            [meta, vcf, index, hap, legend]
        }

    // Check that all input files have the correct index
    checkFileIndex(ch_input.mix(ch_input_truth, ch_ref_gen, ch_panel))

    // Chunk model
    chunk_model = params.chunk_model

    emit:
    input                = ch_input         // [ [meta], file, index ]
    input_truth          = ch_input_truth   // [ [meta], file, index ]
    fasta                = ch_ref_gen       // [ [genome], fasta, fai ]
    panel                = ch_panel         // [ [panel, chr], vcf, index ]
    depth                = ch_depth         // [ [depth], depth ]
    regions              = ch_regions       // [ [chr, region], region ]
    gmap                 = ch_map           // [ [map], map ]
    posfile              = ch_posfile       // [ [panel, chr], vcf, index, hap, legend ]
    chunks               = ch_chunks        // [ [chr], txt ]
    chunk_model          = chunk_model
    versions             = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML

    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    // Check that only genome or fasta is provided
    assert (params.genome == null || params.fasta == null) : "Either --genome or --fasta must be provided"
    assert !(params.genome == null && params.fasta == null) : "Only one of --genome or --fasta must be provided"

    // Check that a steps is provided
    assert params.steps : "A step must be provided"

    // Check that at least one tool is provided
    if (params.steps.split(',').contains("impute")) {
        assert params.tools : "No tools provided"
    }

    // Check that input is provided for all steps, except panelprep
    if (params.steps.split(',').contains("all") || params.steps.split(',').contains("impute") || params.steps.split(',').contains("simulate") || params.steps.split(',').contains("validate")) {
        assert params.input : "No input provided"
    }

    // Check that posfile and chunks are provided when running impute only. Steps with panelprep generate those files.
    if (params.steps.split(',').contains("impute") && !params.steps.split(',').find { it in ["all", "panelprep"] }) {
        // Required by all tools except glimpse2
        if (!params.tools.split(',').find { it in ["glimpse2"] }) {
                assert params.posfile : "No --posfile provided for --steps impute"
        }
        // Required by all tools except STITCH
        if (params.tools != "stitch") {
                assert params.chunks : "No --chunks provided for --steps impute"
        }
        // Required by GLIMPSE1 and GLIMPSE2 only
        if (params.tools.split(',').contains("glimpse")) {
                assert params.panel : "No --panel provided for imputation with GLIMPSE"
        }

        // Check that input_truth is provided when running validate
        if (params.steps.split(',').find { it in ["all", "validate"] } ) {
            assert params.input_truth : "No --input_truth was provided for --steps validate"
        }
    }

    // Emit a warning if both panel and (chunks || posfile) are used as input
    if (params.panel && params.chunks && params.steps.split(',').find { it in ["all", "panelprep"]} ) {
        log.warn("Both `--chunks` and `--panel` have been provided. Provided `--chunks` will override `--panel` generated chunks in `--steps impute` mode.")
    }
    if (params.panel && params.posfile && params.steps.split(',').find { it in ["all", "panelprep"]} ) {
        log.warn("Both `--posfile` and `--panel` have been provided. Provided `--posfile` will override `--panel` generated posfile in `--steps impute` mode.")
    }

    // Emit an info message when using external panel and impute only
    if (params.panel && params.steps.split(',').find { it in ["impute"] } && !params.steps.split(',').find { it in ["all", "panelprep"] } ) {
        log.info("Provided `--panel` will be used in `--steps impute`. Make sure it has been previously prepared with `--steps panelprep`")
    }

    // Emit an error if normalizing step is ignored but samples need to be removed from reference panel
    if (params.steps.split(',').find { it in ["all", "panelprep"] } && params.remove_samples) {
        if (!params.normalize) {
            error("To use `--remove_samples` you need to include `--normalize`.")
        }
    }

    // Check that the chunk model is provided
    assert params.chunk_model : "No chunk model provided"

    return null
}

//
// Check compatibility between input files size, extension and tools
//
def validateInputBatchTools(ch_input, batch_size, extension, tools) {
    ch_input
        .count()
        .map{ nb_input ->
            if (extension ==~ "(vcf|bcf)(.gz)?") {
                if (tools.contains("stitch") || tools.contains("quilt")) {
                    error "Stitch or Quilt software cannot run with VCF or BCF files. Please provide alignment files (i.e. BAM or CRAM)."
                }
                if (nb_input > 1) {
                    error "When using a Variant Calling Format file as input, only one file can be provided. If you have multiple single-sample VCF files, please merge them into a single multisample VCF file."
                }
            }

            if (nb_input > batch_size) {
                if (tools.contains("glimpse2") || tools.contains("quilt")) {
                    log.warn("Glimpse2 or Quilt software is selected and the number of input files (${nb_input}) is less than the batch size (${batch_size}). The input files will be processed in ${Math.ceil(nb_input / batch_size) as int} batches.")
                }
                if (tools.contains("stitch") || tools.contains("glimpse1")) {
                    error "Stitch or Glimpse1 software is selected and the number of input files (${nb_input}) is less than the batch size (${batch_size}). Splitting the input files in batches would induce batch effect."
                }
            }
        }
    return null
}

//
// Check if posfile is compatible with tools and steps selected
//
def validatePosfileTools(ch_posfile, tools, steps){
    ch_posfile
        .map{ _meta, vcf, index, hap, legend ->
            if (tools.contains("glimpse1")) {
                assert legend : "Glimpse1 tool needs a legend file provided in the posfile. This file can be created through the panelprep step."
            }
            if (tools.contains("stitch")) {
                assert legend : "Stitch tool needs a legend file provided in the posfile. This file can be created through the panelprep step."
            }
            if (tools.contains("quilt")) {
                assert legend : "Quilt tool needs a legend file provided in the posfile. This file can be created through the panelprep step."
                assert hap : "Quilt tool needs a hap file provided in the posfile. This file can be created through the panelprep step."
            }
            if (steps.contains("validate")) {
                assert vcf : "Validation step needs a vcf file provided in the posfile for the allele frequency. This file can be created through the panelprep step."
                assert index : "Validation step needs an index file provided in the posfile for the allele frequency. This file can be created through the panelprep step."
            }
        }
    return null
}

//
// Extract contig names from channel meta map
//
def extractChr(ch_input) {
    ch_input.map { [it[0].chr] }
        .collect()
        .toList()
}

//
// Check if all contigs in a are present in b
// Give back the intersection of a and b
//
def checkMetaChr(chr_a, chr_b, name){
    def intersect = chr_a
        .combine(chr_b)
        .map{
            a, b ->
            if (b != [[]] && !(a - b).isEmpty()) {
                def chr_names = (a - b).size() > params.max_chr_names ? (a - b)[0..params.max_chr_names - 1] + ['...'] : (a - b)
                def verb = (a - b).size() == 1 ? "is" : "are"
                log.warn "Chr : ${chr_names} ${verb} missing from ${name}"
                return (a-b)
            }
            return []
        }
        .flatten()
    return intersect
}

//
// Get file extension
//
def getFileExtension(file) {
    def file_name = ""

    if (file instanceof Path || file instanceof nextflow.file.http.XPath) {
        file_name = file.name
    } else if (file instanceof CharSequence) {
        file_name = file.toString()
    } else if (file instanceof List) {
        return file.collect { getFileExtension(it) }
    } else {
        error "Type not supported: ${file.getClass()}"
    }

    // Remove .gz if present and get the last part after splitting by "."
    return file_name.replace(".gz", "").split("\\.").last()
}

//
// Check if all input files have the same extension
//
def getFilesSameExt(ch_input) {
    return ch_input
        .map { getFileExtension(it[1]) } // Extract files extensions
        .toList()  // Collect extensions into a list
        .map { extensions ->
            if (extensions.unique().size() > 1) {
                error "All input files must have the same extension: ${extensions.unique()}"
            }
            return extensions[0]
        }
}

//
// Check correspondance file / index
//
def checkFileIndex(ch_input) {
    ch_input.map {
        meta, file, index ->
        def file_ext = getFileExtension(file)
        def index_ext = getFileExtension(index)
        if (file_ext in ["vcf", "bcf"] &&  !(index_ext in ["tbi", "csi"]) ) {
            log.info("File: ${file} ${file_ext}, Index: ${index} ${index_ext}")
            error "${meta}: Index file for [.vcf, .vcf.gz, bcf] must have the extension [.tbi, .csi]"
        }
        if (file_ext == "bam" && index_ext != "bai") {
            error "${meta}: Index file for .bam must have the extension .bai"
        }
        if (file_ext == "cram" && index_ext != "crai") {
            error "${meta}: Index file for .cram must have the extension .crai"
        }
        if (file_ext in ["fa", "fasta"] && index_ext != "fai") {
            error "${meta}: Index file for [fa, fasta] must have the extension .fai"
        }
    }
    return null
}

//
// Export a channel to a CSV file with correct paths
//
def exportCsv(ch_files, metas, header, name, outdir) {
    ch_files.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/${outdir}") { it ->
        def meta = ""
        def file = ""
        metas.each { i ->
            meta += "${it[0][i]},"
        }
        it[1].each { i ->
            file += "${params.outdir}/${i.value}/${it[i.key].fileName},"
        }
        file = file.substring(0, file.length() - 1) // remove last comma
        ["${name}", "${header}\n${meta}${file}\n"]
    }
    return null
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (meta, bam, bai) = input
    // Check that individual IDs are unique
    // No check for the moment

    return [meta, bam, bai]
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "BCFtools (Danecek et al. 2021),",
        params.tools ? params.tools.split(',').contains("glimpse")   ? "GLIMPSE (Rubinacci et al. 2020)," : "" : "",
        params.tools ? params.tools.split(',').contains("glimpse2")  ? "GLIMPSE2 (Rubinacci et al. 2023)," : "": "",
        params.tools ? params.tools.split(',').contains("quilt")     ? "QUILT (Davies et al. 2021)," : "": "",
        "SAMtools (Li et al. 2009),",
        params.tools ? params.phase ? "SHAPEIT5 (Hofmeister et al. 2023)," : "": "",
        params.tools ? params.phase ? "BEDtools (Quinlan and Hall 2010)," : "": "",
        params.tools ? params.tools.split(',').contains("stitch")    ? "STITCH (Davies et al. 2016)," : "": "",
        "Tabix (Li et al. 2011),",
        params.tools ? params.compute_freq  ? "VCFlib (Garrison et al. 2022)," : "": "",
        "."
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        params.phase ? "<li>Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi:10.1093/bioinformatics/btq033.</li>": "",
        "<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi:10.1093/bioinformatics/btp352.</li>",
        "<li>Li H. (2011). Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics. 2011 Mar 1;27(5):718-9. doi:10.1093/bioinformatics/btq671.</li>",
        params.tools ? params.tools.split(',').contains("quilt") ? "<li>Davies RW, Kucka M, Su D, Shi S, Flanagan M, Cunniff CM, Chan YF, & Myers S. (2021). Rapid genotype imputation from sequence with reference panels. Nature Genetics. doi:10.1038/s41588-021-00877-0.</li>" : "": "",
        params.tools ? params.tools.split(',').contains("glimpse") ? "<li>Rubinacci S, Ribeiro DM, Hofmeister RJ, & Delaneau O. (2021). Efficient phasing and imputation of low-coverage sequencing data using large reference panels. Nature Genetics. doi:10.1038/s41588-020-00756-0.</li>" : "": "",
        params.tools ? params.tools.split(',').contains("glimpse2") ? "<li>Rubinacci S, Hofmeister RJ, Sousa da Mota B, & Delaneau O. (2023). Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. Nature Genetics. doi:10.1038/s41588-023-01438-3.</li>" : "": "",
        params.phase ? "<li>Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O. (2023). Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. Nat Genet. 2023 Jul;55(7):1243-1249. doi:10.1038/s41588-023-01415-w.</li>" : "",
        params.tools ? params.tools.split(',').contains("stitch") ? "<li>Davies RW, Flint J, Myers S, & Mott R. (2016). Rapid genotype imputation from sequence without reference panels. Nature Genetics.</li>" : "": "",
        params.compute_freq ? "<li>Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P. (2022). A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. PLoS Comput Biol 18(5).</li>" : "",
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

