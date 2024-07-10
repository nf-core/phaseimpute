//
// Subworkflow with functionality specific to the nf-core/phaseimpute pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'
include { GET_REGION                } from '../get_region'
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

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
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
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
        .fromSamplesheet("input")
        .map {
            meta, file, index ->
                [ meta, file, index ]
        }
        // Check if all extension are identical
        getAllFilesExtension(ch_input)
    } else {
        ch_input = Channel.of([[], [], []])
    }
    //
    // Create channel from input file provided through params.input_truth
    //
    if (params.input_truth) {
        if (params.input_truth.endsWith("csv")) {
            ch_input_truth = Channel
                .fromSamplesheet("input_truth")
                .map {
                    meta, file, index ->
                        [ meta, file, index ]
                }
            // Check if all extension are identical
            getAllFilesExtension(ch_input_truth)
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
            print("Panel file provided as input is a samplesheet")
            ch_panel = Channel.fromSamplesheet("panel")
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
        ch_regions = Channel.fromSamplesheet("input_region")
            .map{ chr, start, end -> [["chr": chr], chr + ":" + start + "-" + end]}
            .map{ metaC, region -> [metaC + ["region": region], region]}
    } else {
        error "Region file provided is of another format than CSV (not yet supported). Please separate your reference genome by chromosome and use the samplesheet format."
    }

    //
    // Create map channel
    //
    if (params.map) {
        if (params.map.endsWith(".csv")) {
            print("Map file provided as input is a samplesheet")
            ch_map = Channel.fromSamplesheet("map")
        } else {
            error "Map file provided is of another format than CSV (not yet supported). Please separate your reference genome by chromosome and use the samplesheet format."
        }
    } else {
        ch_map = ch_regions
            .map{ metaCR, regions -> [metaCR.subMap("chr"), []] }
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
        ch_posfile = Channel
                .fromSamplesheet("posfile") // ["panel", "chr", "vcf", "index", "hap", "legend"]
    } else {
        ch_posfile = Channel.of([[],[],[],[],[]])
    }

    //
    // Create chunks channel
    //
    if (params.chunks) {
        ch_chunks = Channel
            .fromSamplesheet("chunks")
    } else {
        ch_chunks = Channel.of([[],[]])
    }

    //
    // Check contigs name in different meta map
    //
    // Collect all chromosomes names in all different inputs
    chr_ref = ch_ref_gen.map { meta, fasta, fai -> [fai.readLines()*.split('\t').collect{it[0]}] }
    chr_regions = extractChr(ch_regions)

    // Check that the chromosomes names that will be used are all present in different inputs
    checkChr(chr_regions, chr_ref, "reference genome")
    checkChr(chr_regions, extractChr(ch_chunks), "chromosome chunks")
    checkChr(chr_regions, extractChr(ch_map), "genetic map")
    checkChr(chr_regions, extractChr(ch_panel), "reference panel")
    checkChr(chr_regions, extractChr(ch_posfile), "position")

    // Check that all input files have the correct index
    checkFileIndex(ch_input)
    checkFileIndex(ch_input_truth)
    checkFileIndex(ch_ref_gen)
    checkFileIndex(ch_panel)

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
    versions             = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
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

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
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
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    // Check that only genome or fasta is provided
    assert params.genome == null || params.fasta == null, "Either --genome or --fasta must be provided"
    assert !(params.genome == null && params.fasta == null), "Only one of --genome or --fasta must be provided"

    // Check that a steps is provided
    assert params.steps, "A step must be provided"

    // Check that at least one tool is provided
    if (params.steps.split(',').contains("impute")) {
        assert params.tools, "No tools provided"
    }

    // Check that input is provided for all steps, except panelprep
    if (params.steps.split(',').contains("all") || params.steps.split(',').contains("impute") || params.steps.split(',').contains("simulate") || params.steps.split(',').contains("validate")) {
        assert params.input, "No input provided"
    }

    // Check that posfile and chunks are provided when running impute only. Steps with panelprep generate those files.
    if (params.steps.split(',').contains("impute") && !params.steps.split(',').find { it in ["all", "panelprep"] }) {
        // Required by all tools except glimpse2
        if (!params.tools.split(',').find { it in ["glimpse2"] }) {
                assert params.posfile, "No --posfile provided for --steps impute"
        }
        // Required by all tools except STITCH
        if (params.tools != "stitch") {
                assert params.chunks, "No --chunks provided for --steps impute"
        }
        // Required by GLIMPSE1 and GLIMPSE2 only
        if (params.tools.split(',').contains("glimpse")) {
                assert params.panel, "No --panel provided for imputation with GLIMPSE"
        }

        // Check that input_truth is provided when running validate
        if (params.steps.split(',').find { it in ["all", "validate"] } ) {
            assert params.input_truth, "No --input_truth was provided for --steps validate"
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
}

//
// Extract contig names from channel meta map
//
def extractChr(channel) {
    channel.map { [it[0].chr] }
        .collect()
        .toList()
}

//
// Check if all contigs in a are present in b
//
def checkChr(chr_a, chr_b, name){
    chr_a
        .combine(chr_b)
        .map{
            a, b ->
            if (b != [[]] && !(a - b).isEmpty()) {
                error "Chr : ${a - b} is missing from ${name}"
            }
        }
}

//
// Get file extension
//
def getFileExtension(file) {
    if (file instanceof String) {
        return file.replace(".gz","").split("\\.").last()
    } else if (file instanceof Path) {
        return file.getName().replace(".gz","").split("\\.").last()
    } else if (file instanceof ArrayList) {
        if (file == []) {
            return null
        } else {
            error "Array not supported"
        }
    } else {
        error "Type not supported: ${file.getClass()}"
    }
}

//
// Check if all input files have the same extension
//
def getAllFilesExtension(ch_input) {
    files_ext = ch_input
        .map { getFileExtension(it[1]) } // Extract files extensions
        .toList()  // Collect extensions into a list
        .map { extensions ->
            if (extensions.unique().size() != 1) {
                error "All input files must have the same extension: ${extensions.unique()}"
            }
            return extensions[0]
        }
}

//
// Check correspondance file / index
//
def checkFileIndex(ch_input) {
    ch_input
        .subscribe {
            meta, file, index ->
            file_ext = getFileExtension(file)
            index_ext = getFileExtension(index)
            log.info("${meta} ${file}, ${index}, ${file_ext}, ${index_ext}")
            if (file_ext in ["vcf", "bcf"] &&  !(index_ext in ["tbi", "csi"]) ) {
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
        meta = ""
        file = ""
        for (i in metas) {
            meta += "${it[0][i]},"
        }
        for (i in it[1]) {
            file += "${params.outdir}/${i.value}/${it[i.key].fileName},"
        }
        file=file.substring(0, file.length() - 1) // remove last comma
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
    // no validation for the moment
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
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
