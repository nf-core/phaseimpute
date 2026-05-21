//
// Subworkflow with functionality specific to the nf-core/phaseimpute pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    ch_ref_gen        // channel: Channel of reference genome file [meta, fasta, fai, gzi]
    sheets_given      //     map: Input sheets given
    steps             //   array: Steps to perform
    tools             //   array: Tools to be used
    max_chr_names     // integer: Maximum number of chromosome to log in warnings and errors
    params_simulate   //     map: Parameters use for simulation step [depth: float, genotype: path]
    params_panelprep  //     map: Parameters use for panelprep step  [normalize: boolean, remove_samples: string, compute_freq: boolean, phase: boolean, chunk_model: string ]
    params_impute     //     map: Parameters use for imputation step [batch_size: integer, k_val: integer, n_gen: integer, buffer: integer]

    main:

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

    def after_text = ""
    def before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/phaseimpute ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/phaseimpute/blob/main/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )


    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    if (tools) {
        log.info("\033[1;37mExtra informations\033[0m")
        log.info("\033[0;34m  Tools selected to be run  :\033[0;32m " + tools.join(",") + "\033[0m")
        log.info("-\033[2m----------------------------------------------------\033[0m-")
    }

    // Initialize variable

    def sheet_target  = sheets_given["input_target"]
    def sheet_truth   = sheets_given["input_truth"]
    def sheet_panel   = sheets_given["input_panel"]
    def sheet_posfile = sheets_given["input_posfile"]
    def sheet_chunks  = sheets_given["input_chunks"]
    def sheet_map     = sheets_given["input_map"]
    def sheet_region  = sheets_given["input_region"]

    def depth = params_simulate["depth"]
    def genotype = params_simulate["genotype"]

    def normalize = params_panelprep["normalize"]
    def remove_samples = params_panelprep["remove_samples"]
    def chunk_model = params_panelprep["chunk_model"]

    def batch_size = params_impute["batch_size"]

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(
        sheet_target, sheet_truth, sheet_panel,
        sheet_posfile, sheet_chunks,
        genotype, remove_samples, normalize,
        chunk_model,
        steps, tools
    )

    //
    // Create channel from input file provided through input
    //
    if (sheet_target) {
        ch_input_target = channel
            .fromList(samplesheetToList(sheet_target, "${projectDir}/assets/schema_input.json"))
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map { meta, file, index ->
                def new_meta = meta + [id:meta.id.toString()]
                [ new_meta + [batch: 0], file, index ]
            } // Set batch to 0 by default
    } else {
        ch_input_target = channel.of([[], [], []])
    }

    // Check that the batch size and extension is compatible with the tools
    validateInputBatchTools(
        ch_input_target,
        batch_size,
        getFilesSameExt(ch_input_target),
        tools
    )

    //
    // Create channel from input file provided through input_truth
    //
    if (sheet_truth) {
        if (sheet_truth.endsWith("csv")) {
            ch_input_truth = channel
                .fromList(samplesheetToList(sheet_truth, "${projectDir}/assets/schema_input.json"))
                .map {
                    meta, file, index ->
                        [ meta + [id:meta.id.toString()], file, index ]
                }
            // Check if all extension are identical
            input_truth_ext = getFilesSameExt(ch_input_truth)
        } else {
            // #TODO Wait for `oneOf()` to be supported in the nextflow_schema.json
            error "Panel file provided is of another format than CSV (not yet supported). Please separate your panel by chromosome and use the samplesheet format."
        }
    } else {
        ch_input_truth = channel.of([[], [], []])
        input_truth_ext = ""
    }

    //
    // Create channel from region input
    //
    if (sheet_region == null){
        // #TODO Add support for string input
        ch_regions = getRegionFromFai("all", ch_ref_gen)
    }  else  if (sheet_region.endsWith(".csv")) {
        log.info "Region file provided as input is a samplesheet"
        ch_regions = channel.from(samplesheetToList(
            sheet_region, "${projectDir}/assets/schema_input_region.json"
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
    // Create channel for panel
    //
    if (sheet_panel) {
        if (sheet_panel.endsWith("csv")) {
            log.info "Panel file provided as input is a samplesheet"
            ch_panel = channel.fromList(samplesheetToList(
                sheet_panel, "${projectDir}/assets/schema_input_panel.json"
            )).map {
                meta, file, index ->
                    [ meta + [panel_id:meta.panel_id.toString()], file, index ]
            }
        } else {
            // #TODO Wait for `oneOf()` to be supported in the nextflow_schema.json
            error "Panel file provided is of another format than CSV (not yet supported). Please separate your panel by chromosome and use the samplesheet format."
        }
    } else {
        ch_panel = ch_regions
            .map{ metaCR, _regions -> [[panel_id: "None"] + metaCR.subMap("chr"), [], []] }
    }

    //
    // Create map channel
    //
    if (sheet_map) {
        if (sheet_map.endsWith(".csv")) {
            log.info "Map file provided as input is a samplesheet"
            ch_map = channel.fromList(samplesheetToList(sheet_map, "${projectDir}/assets/schema_map.json"))
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
    if (depth) {
        ch_depth = channel.of([[depth: depth], depth])
    } else {
        ch_depth = channel.of([[],[]])
    }

    //
    // Create genotype array channel
    //
    if (genotype) {
        ch_genotype = channel.of([[gparray: genotype], genotype])
    } else {
        ch_genotype = channel.of([[],[]])
    }

    //
    // Create posfile channel
    //
    if (sheet_posfile) {
        ch_posfile = channel // ["meta", "vcf", "index", "hap", "legend", "posfile"]
            .fromList(samplesheetToList(sheet_posfile, "${projectDir}/assets/schema_posfile.json"))
            .map { meta, vcf, index, hap, legend, posfile ->
                [ meta + [panel_id:meta.panel_id.toString()], vcf, index, hap, legend, posfile ]
            }
    } else {
        ch_posfile = ch_panel
            .map{ metaPC, _vcf, _index -> [metaPC, [], [], [], [], []]}
    }

    if (!steps.contains("panelprep") & !steps.contains("all")) {
        validatePosfileTools(
            ch_posfile,
            tools,
            steps,
            input_truth_ext
        )
    }

    //
    // Create chunks channel
    //
    if (sheet_chunks) {
        ch_chunks = channel
            .fromList(samplesheetToList(sheet_chunks, "${projectDir}/assets/schema_chunks.json"))
            .map { meta, chunks ->
                [ meta + [panel_id:meta.panel_id.toString()], chunks ]
            }
    } else {
        ch_chunks = ch_panel
            .map{ metaPC, _vcf, _index -> [metaPC, []] }
    }

    //
    // Check panel, chunks and posfile have same panel id
    //
    panel_panelid   = ch_panel.map{ metaPC, _vcf, _index -> [metaPC.panel_id]}.unique()
    chunks_panelid  = ch_chunks.map{ metaPC, _chunks -> [metaPC.panel_id]}.unique()
    posfile_panelid = ch_posfile.map{ metaPC, _vcf, _index, _hap, _legend, _posfile -> [metaPC.panel_id]}.unique()

    // Get all unique panel id except None
    panel_id = panel_panelid
        .mix(chunks_panelid, posfile_panelid)
        .flatten()
        .filter { it -> it != "None" }
        .unique()

    // Check uniqueness of panel_id
    // TODO add support for multiple panel
    panel_id
        .collect()
        .map{ panel_ids ->
            if (panel_ids.size() != 1) {
                error "Multiple panel IDs detected: ${panel_ids}. Please provide only one across panel, chunks and posfile."
            }
        }

    //
    // Check contigs name in different meta map
    //
    // Collect all chromosomes names in all different inputs
    chr_ref = ch_ref_gen.map {
        _meta, _fasta, fai_file, _gzi_file -> [
            fai_file.readLines()*.split('\t').collect{cols -> cols[0]}
        ]
    }
    chr_regions = extractChr(ch_regions)

    // Check that the chromosomes names that will be used are all present in different inputs
    chr_ref_mis     = checkMetaChr(chr_regions, chr_ref, "reference genome", max_chr_names)
    chr_chunks_mis  = checkMetaChr(chr_regions, extractChr(ch_chunks), "chromosome chunks", max_chr_names)
    chr_map_mis     = checkMetaChr(chr_regions, extractChr(ch_map), "genetic map", max_chr_names)
    chr_panel_mis   = checkMetaChr(chr_regions, extractChr(ch_panel), "reference panel", max_chr_names)
    chr_posfile_mis = checkMetaChr(chr_regions, extractChr(ch_posfile), "position", max_chr_names)

    // Compute the intersection of all chromosomes names
    chr_all_mis_ch = chr_ref_mis.concat(chr_chunks_mis, chr_map_mis, chr_panel_mis, chr_posfile_mis)
        .unique()
    chr_all_mis = chr_all_mis_ch.toList()

    chr_all_mis_for_combine = chr_all_mis.map { chr -> [chr] }

    chr_all_mis.subscribe{ chr ->
        if (chr.size() > 0) {
            def chr_names = chr.size() > max_chr_names ? chr[0..max_chr_names - 1] + ['...'] : chr
            log.warn "The following contigs are absent from at least one file : ${chr_names} and therefore won't be used"
        }
    }

    ch_regions = ch_regions
        .combine(chr_all_mis_for_combine)
        .filter { meta, _regions, chr_mis ->
            !(meta.chr in chr_mis)
        }
        .map { meta, regions, _chr_mis -> [meta, regions] }
        .ifEmpty { error "No regions left to process" }

    ch_regions
        .map { _metaCR, region -> region }
        .collect()
        .subscribe { region -> log.info "The following contigs will be processed: ${region}" }

    // Remove other contigs from panel and posfile files
    ch_panel = ch_panel
        .combine(ch_regions.collect{ metaCR, _region -> metaCR.chr }.toList())
        .filter { meta, _vcf, _index, chrs ->
            meta.chr in chrs
        }
        .map {meta, vcf, index, _chrs ->
            [meta, vcf, index]
        }

    ch_posfile = ch_posfile
        .combine(ch_regions.collect{ metaCR, _region -> metaCR.chr }.toList())
        .filter { meta, _vcf, _index, _hap, _legend, _posfile, chrs ->
            meta.chr in chrs
        }
        .map {meta, vcf, index, hap, legend, posfile, _chrs ->
            [meta, vcf, index, hap, legend, posfile]
        }

    // Combine map and panel for joint operations
    ch_map = ch_map
        .combine(ch_panel.map{ metaPC, _vcf, _index -> [
            metaPC.subMap("chr"), metaPC
        ]}, by: 0)
        .map{ _metaC, map, metaPC -> [
            metaPC, map
        ]}

    // Check that all input files have the correct index
    checkFileIndex(ch_input_target.mix(ch_input_truth, ch_ref_gen, ch_panel))

    // Make available both index if present
    ch_fasta_index = ch_ref_gen.map{ meta, fasta, fai, gzi -> [
        meta, fasta, gzi ? [fai, gzi] : [fai]
    ]}

    emit:
    ch_input_target = ch_input_target // [ [meta], file, index ]
    ch_input_truth  = ch_input_truth  // [ [meta], file, index ]
    ch_fasta_index  = ch_fasta_index  // [ [genome], fasta, [fai, gzi] ]
    ch_panel        = ch_panel        // [ [panel_id, chr], vcf, index ]
    ch_depth        = ch_depth        // [ [depth], depth ]
    ch_regions      = ch_regions      // [ [chr, region], region ]
    ch_map          = ch_map          // [ [chr], map ]
    ch_posfile      = ch_posfile      // [ [panel_id, chr], vcf, index, hap, legend, posfile ]
    ch_chunks       = ch_chunks       // [ [panel_id, chr], txt ]
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

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
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
def validateInputParameters(
    sheet_target, sheet_truth, sheet_panel,
    sheet_posfile, sheet_chunks,
    genotype, remove_samples, normalize,
    chunk_model,
    steps, tools
) {
    // Check that a steps is provided
    if (!steps) {
        error "A step must be provided"
    }

    // Check that genotype simulation is only used with simulate/all steps
    if (genotype && !steps.contains("impute") && !steps.contains("all")) {
        error("`--genotype` is only supported with `--steps simulate` or `--steps all`. Genotype simulation is not yet implemented.")
    }

    // Check that at least one tool is provided
    if (steps.contains("impute")) {
        if (!tools) {
            error "No tools provided"
        }
    }

    // Check that input is provided for all steps, except panelprep
    if (steps.contains("all") || steps.contains("impute") || steps.contains("simulate") || steps.contains("validate")) {
        if (!sheet_target) {
            error "No input provided"
        }
    }

    // Check that posfile and panel are provided when running impute only
    if (steps.contains("impute") && !steps.find { step -> step in ["all", "panelprep"] }) {
        // Required by all tools except glimpse2, quilt2, beagle5, minimac4
        if (!tools.find { tool -> tool in ["glimpse2", "quilt2", "beagle5", "minimac4"] }) {
            if (!sheet_posfile) {
                error "No --posfile provided for --steps impute"
            }
        }
        // Required by panel-backed imputation tools
        if (tools.find { tool -> tool in ["glimpse1", "glimpse2", "quilt2"] }) {
            if (!sheet_panel) {
                error "No --panel provided for imputation with GLIMPSE1, GLIMPSE2 or QUILT2"
            }
        }
    }

    // Check that input_truth is provided when running validate
    if (steps.find { step -> step in ["validate"] } && !steps.find { step -> step in ["simulate"] }) {
        if (!sheet_truth) {
            error "No --input_truth was provided for --steps validate"
        }
    }

    // Emit a warning if both panel and (chunks || posfile) are used as input
    if (sheet_panel && sheet_chunks && steps.find { step -> step in ["all", "panelprep"]} ) {
        log.warn("Both `--chunks` and `--panel` have been provided. Provided `--chunks` will override `--panel` generated chunks in `--steps impute` mode.")
    }
    if (sheet_panel && sheet_posfile && steps.find { step -> step in ["all", "panelprep"]} ) {
        log.warn("Both `--posfile` and `--panel` have been provided. Provided `--posfile` will override `--panel` generated posfile in `--steps impute` mode.")
    }

    // Emit an info message when using external panel and impute only
    if (sheet_panel && steps.find { step -> step in ["impute"] } && !steps.find { step -> step in ["all", "panelprep"] } ) {
        log.info("Provided `--panel` will be used in `--steps impute`. Make sure it has been previously prepared with `--steps panelprep`")
    }

    // Emit an error if normalizing step is ignored but samples need to be removed from reference panel
    if (steps.find { step -> step in ["all", "panelprep"] } && remove_samples) {
        if (!normalize) {
            error("To use `--remove_samples` you need to include `--normalize`.")
        }
    }

    // Check that the chunk model is provided
    if (!chunk_model) {
        error "No chunk model provided"
    }

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
                if (tools.contains("stitch") || tools.contains("quilt") || tools.contains("quilt2")) {
                    error "Stitch, QUILT and QUILT2 software cannot run with VCF or BCF files. Please provide alignment files (i.e. BAM or CRAM)."
                }
                if (nb_input > 1) {
                    error "When using a Variant Calling Format file as input, only one file can be provided. If you have multiple single-sample VCF files, please merge them into a single multisample VCF file."
                }
            }

            if (extension ==~ "(bam|cram)?") {
                if (tools.contains("beagle5") || tools.contains("minimac4")) {
                    error "Beagle5 and Minimac4 softwares cannot run with BAM or CRAM alignement files. Please provide variant calling format files (i.e. VCF or BCF)."
                }
            }

            if (nb_input > batch_size) {
                if (tools.contains("glimpse2") || tools.contains("quilt") || tools.contains("quilt2")) {
                    log.warn("Glimpse2, QUILT or QUILT2 software is selected and the number of input files (${nb_input}) is less than the batch size (${batch_size}). The input files will be processed in ${Math.ceil(nb_input / batch_size) as int} batches.")
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
def validatePosfileTools(ch_posfile, tools, steps, truth_extension){
    ch_posfile
        .map{ _meta, vcf, index, hap, legend, posfile ->
            if (tools.contains("glimpse1")) {
                assert posfile : "Glimpse1 tool needs a posfile file with CHROM\tPOS\tREF,ALT columns. This file can be created through the panelprep step."
            }
            if (tools.contains("stitch")) {
                assert posfile : "You have not provided a posfile and you've requested to use STITCH. In this pipeline, using STITCH requires a posfile file with CHROM\tPOS\tREF,ALT columns. This file is generated automatically in the panelprep step."
            }
            if (tools.contains("quilt")) {
                assert legend : "Quilt tool needs a legend file provided in the posfile. This file can be created through the panelprep step."
                assert hap : "Quilt tool needs a hap file provided in the posfile. This file can be created through the panelprep step."
            }
            if (steps.contains("validate")) {
                assert vcf : "Validation step needs a vcf file provided in the posfile for the allele frequency. This file can be created through the panelprep step."
                assert index : "Validation step needs an index file provided in the posfile for the allele frequency. This file can be created through the panelprep step."
                if (truth_extension =~ "bam|cram"){
                    assert posfile : "You have not provided a posfile and you've requested to use the validation step with bam files. This step requires a posfile file with CHROM\tPOS\tREF,ALT columns to call the variants from the truth BAM file. This file is generated automatically in the panelprep step."
                }
            }
        }

    ch_posfile
        .map{ _meta, _vcf, _index, _hap, _legend, posfile ->
            if (posfile) {
                def lines = []
                def pathFile = posfile instanceof String ? file(posfile) : posfile
                pathFile.withInputStream { stream ->
                    def reader = pathFile.name.endsWith('.gz') ?
                        new java.util.zip.GZIPInputStream(stream).newReader() :
                        new InputStreamReader(stream)

                    reader.withReader { r ->
                        (1..3).each {
                            def line = r.readLine()
                            if (line != null) lines << line
                        }
                    }
                }

                // Validate first 3 lines (or fewer if file is shorter)
                lines.each { line ->
                    def fields = line.split("\t")
                    assert fields.size() == 3 : "Expected 3 columns in ${posfile.name}, found ${fields.size()} in line: ${line}"
                    assert fields[2].contains(",") : "Third column must contain comma in ${posfile.name}, line: ${line}"
                }
            }
        }

    return null
}

//
// Extract contig names from channel meta map
//
def extractChr(ch_input) {
    ch_input.map { it -> [it[0].chr] }
        .collect()
        .toList()
}

//
// Check if all contigs in a are present in b
// Give back the intersection of a and b
//
def checkMetaChr(chr_a, chr_b, name, max_chr_names){
    def intersect = chr_a
        .combine(chr_b)
        .map{
            a, b ->
            if (b != [[]] && !(a - b).isEmpty()) {
                def chr_names = (a - b).size() > max_chr_names ? (a - b)[0..max_chr_names - 1] + ['...'] : (a - b)
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
// Get region from fasta fai file
//
def getRegionFromFai(region_selected, ch_fasta) {
    def ch_regions = channel.empty()
    // Gather regions to use and create the meta map
    if (region_selected ==~ '^(chr)?[0-9XYM]+$' || region_selected == "all") {
        ch_regions = ch_fasta.map{ _meta, _fasta, fai, _gzi -> fai}
            .splitCsv(header: ["chr", "size", "offset", "lidebase", "linewidth", "qualoffset"], sep: "\t")
            .map{it -> [chr:it.chr, region:"0-"+it.size]}
        if (region_selected != "all") {
            ch_regions = ch_regions.filter{ it -> it.chr == region_selected}
        }
        ch_regions = ch_regions
            .map{ it -> [[chr: it.chr, region: it.chr + ":" + it.region], it.chr + ":" + it.region]}
    } else {
        if (region_selected ==~ '^chr[0-9XYM]+:[0-9]+-[0-9]+$') {
            ch_regions = channel.from([region_selected])
                .map{ it -> [[chr: it.split(":")[0], "region": it], it]}
        } else {
            error "Invalid region_selected: ${region_selected}"
        }
    }
    return ch_regions
}

//
// Get file extension
//
def getFileExtension(file) {
    def file_name = ""
    if (file instanceof Path) {
        file_name = file.name
    } else if (file instanceof java.net.URL) {
        file_name = file.path.tokenize('/')[-1]
    } else if (file instanceof CharSequence) {
        file_name = file.toString()
    } else if (file instanceof List) {
        return file.collect { it -> getFileExtension(it) }
    } else if (file == null) {
        return ""
    } else {
        error "Type not supported: ${file.getClass()}"
    }
    // Remove .gz if present at the end and get the last part after splitting by "."
    return file_name.replaceAll(/\.gz$/, "").tokenize('.').last()
}

//
// Check if all input files have the same extension
//
def getFilesSameExt(ch_input) {
    return ch_input
        .map { it -> getFileExtension(it[1]) } // Extract files extensions
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
        tuple ->
        def meta = tuple[0]
        def file = tuple[1]
        def index = tuple[2]
        def gzi = tuple.size() > 3 ? tuple[3] : null

        def file_ext = getFileExtension(file)
        def index_ext = getFileExtension(index)
        if (file_ext in ["vcf", "bcf"] &&  !(index_ext in ["tbi", "csi"]) ) {
            log.info("File: ${file} ${file_ext}, Index: ${index} ${index_ext}")
            error "${meta}: Index file for .vcf, .vcf.gz and .bcf must have the extension .tbi or .csi"
        }
        if (file_ext == "bam" && !(index_ext in ["bai", "csi"])) {
            error "${meta}: Index file for .bam must have the extension .bai or .csi"
        }
        if (file_ext == "cram" && index_ext != "crai") {
            error "${meta}: Index file for .cram must have the extension .crai"
        }
        if (file_ext in ["fa", "fasta"] && index_ext != "fai") {
            error "${meta}: Index file for .fa and .fasta must have the extension .fai"
        }

        def file_name = file.toString()
        if (file_ext in ["fa", "fasta"] && (file_name.endsWith('.gz') || file_name.endsWith('.bgz'))) {
            if (getFileExtension(gzi) != "gzi") {
                error "${meta}: Index file for .(fa,fasta).(gz,bgz) must have the extension .gzi"
            }
        }
    }
    return null
}

//
// Export a channel to a CSV file with correct paths
//
def exportCsv(ch_files, metas, header, name, outdir, folder) {
    ch_files.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/${folder}") { it ->
        def meta = ""
        def file = ""
        metas.each { i ->
            meta += "${it[0][i]},"
        }
        it[1].each { i ->
            file += "${outdir}/${i.value}/${it[i.key].fileName},"
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
// Generate methods description for MultiQC
//
def toolCitationText(steps, tools, normalize, remove_samples, compute_freq, phase) {
    def tool_citation = [
        BEAGLE5 : "Beagle5 (Browning et al. 2018)",
        BCFTOOLS: "BCFtools (Danecek et al. 2021)",
        SAMTOOLS: "SAMtools (Danecek et al. 2021)",
        MINIMAC4: "Minimac4 (Das et al. 2016)",
        STITCH  : "STITCH (Davies et al. 2016)",
        QUILT   : "QUILT (Davies et al. 2021)",
        QUILT2  : "QUILT2 (Li et al. 2026)",
        MULTIQC : "MultiQC (Ewels et al. 2016)",
        VCFLIB  : "vcflib (Garrison et al. 2022)",
        SHAPEIT5: "SHAPEIT5 (Hofmeister et al. 2023)",
        TABIX   : "Tabix (Li H et al. 2011)",
        GLIMPSE1: "GLIMPSE (Rubinacci et al. 2021)",
        GLIMPSE2: "GLIMPSE2 (Rubinacci et al. 2023)",
    ]

    if (steps.contains("all")) {
        steps = ["simulate", "panelprep", "impute", "validate"]
    }

    def text_simulate = [
        "Low-coverage sequencing data simulation was performed with",
        "${tool_citation.SAMTOOLS} subcommand 'depth' and 'view' for downsampling high-coverage BAM files."
    ].join(' ').trim()

    def text_panelprep = [
        "Reference panel preparation followed several steps.",
        normalize && remove_samples ? "The reference panel genotypes were normalized and samples: " + remove_samples.split(",").join(", ") + " were removed" :
            normalize ? "The reference panel genotypes were normalized" :
                remove_samples ? "Samples: " + remove_samples.split(",").join(", ") + " were removed from the reference panel genotypes" :
                    "No normalization or sample removal were performed on the reference panel genotypes.",
        normalize || remove_samples ? "followed by site extraction and format conversion using ${tool_citation.BCFTOOLS}.":
            "Site extraction and format conversion was done using ${tool_citation.BCFTOOLS}.",
        compute_freq ? "Allele frequencies were then computed with ${tool_citation.VCFLIB}." : "",
        phase ? "Genotype phasing was performed with ${tool_citation.SHAPEIT5}." : "",
        "Finally, the reference panel was split into per-chromosome chunks using ${tool_citation.GLIMPSE1}",
        "and ${tool_citation.GLIMPSE2}."
    ].join(' ').trim()

    def text_impute = [
        tools.size() > 0 ? tools.size() == 1 ? "Imputation tool used was:" :
            "Imputation tools used were:" : "",
        [
            tools.contains("glimpse1")    ? "${tool_citation.GLIMPSE1}" +
                " with variants called using ${tool_citation.BCFTOOLS} mpileup followed by indexation with ${tool_citation.TABIX}" +
                " when BAM files were provided" : "",
            tools.contains("glimpse2")   ? "${tool_citation.GLIMPSE2}" : "",
            tools.contains("quilt")      ? "${tool_citation.QUILT}"    : "",
            tools.contains("quilt2")     ? "${tool_citation.QUILT2}"   : "",
            tools.contains("stitch")     ? "${tool_citation.STITCH}"   : "",
            tools.contains("beagle5")    ? "${tool_citation.BEAGLE5}"  : "",
            tools.contains("minimac4")   ? "${tool_citation.MINIMAC4}" : ""
        ].findAll{ it -> it != "" }.join(', ') + "."
    ].join(' ').trim()

    def text_validate = [
        "Imputation accuracy was assessed by comparing imputed genotypes to truth data using ${tool_citation.GLIMPSE2}.",
        "Truth genotypes were obtained either from array genotyping data provided as input or from high-coverage sequencing data from which",
        "genotypes were called using ${tool_citation.BCFTOOLS} mpileup followed by indexation with ${tool_citation.TABIX}."
    ].join(' ').trim()

    def text_multiqc = "Pipeline results statistics were summarised with ${tool_citation.MULTIQC}."

    def citation_text = [
        "Tools used in the workflow included the following.",
        steps.contains("simulate")  ? text_simulate  : "",
        steps.contains("panelprep") ? text_panelprep : "",
        steps.contains("impute")    ? text_impute    : "",
        steps.contains("validate")  ? text_validate  : "",
        text_multiqc
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText(steps, tools, compute_freq, phase) {
    def tool_biblio = [
        BEAGLE5 : '<li>Browning, B.L., Zhou, Y., Browning, S.R., 2018. A One-Penny Imputed Genome from Next-Generation Reference Panels. Am J Hum Genet 103, 338-348. doi: <a href="https://doi.org/10.1016/j.ajhg.2018.07.015">10.1016/j.ajhg.2018.07.015</a></li>',
        SAM_BCFTOOLS: '<li>Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H., 2021. Twelve years of SAMtools and BCFtools. GigaScience 10, giab008. doi: <a href="https://doi.org/10.1093/gigascience/giab008">10.1093/gigascience/giab008</a></li>',
        MINIMAC4: '<li>Das, S., Forer, L., Schonherr, S., Sidore, C., Locke, A.E., Kwong, A., Vrieze, S.I., Chew, E.Y., Levy, S., McGue, M., Schlessinger, D., Stambolian, D., Loh, P.-R., Iacono, W.G., Swaroop, A., Scott, L.J., Cucca, F., Kronenberg, F., Boehnke, M., Abecasis, G.R., Fuchsberger, C., 2016. Next-generation genotype imputation service and methods. Nat Genet 48, 1284-1287. doi: <a href="https://doi.org/10.1038/ng.3656">10.1038/ng.3656</a></li>',
        STITCH  : '<li>Davies, R.W., Flint, J., Myers, S., Mott, R., 2016. Rapid genotype imputation from sequence without reference panels. Nat Genet 48, 965-969. doi: <a href="https://doi.org/10.1038/ng.3594">10.1038/ng.3594</a></li>',
        QUILT   : '<li>Davies, R.W., Kucka, M., Su, D., Shi, S., Flanagan, M., Cunniff, C.M., Chan, Y.F., Myers, S., 2021. Rapid genotype imputation from sequence with reference panels. Nat Genet 53, 1104-1111. doi: <a href="https://doi.org/10.1038/s41588-021-00877-0">10.1038/s41588-021-00877-0</a></li>',
        QUILT2  : '<li>Li, Z., Albrechtsen, A., Davies, R.W., 2026. Flexible read-aware genotype imputation from sequence using biobank sized reference panels. Nat Commun 17, 524. doi: <a href="https://doi.org/10.1038/s41467-025-67218-1">10.1038/s41467-025-67218-1</a></li>',
        MULTIQC : '<li>Ewels, P., Magnusson, M., Lundin, S., Kaller, M., 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047-3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>',
        VCFLIB  : '<li>Garrison, E., Kronenberg, Z.N., Dawson, E.T., Pedersen, B.S., Prins, P., 2022. A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. PLOS Computational Biology 18, e1009123. doi: <a href="https://doi.org/10.1371/journal.pcbi.1009123">10.1371/journal.pcbi.1009123</a></li>',
        SHAPEIT5: '<li>Hofmeister, R.J., Ribeiro, D.M., Rubinacci, S., Delaneau, O., 2023. Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. Nat Genet 1-7. doi: <a href="https://doi.org/10.1038/s41588-023-01415-w">10.1038/s41588-023-01415-w</a></li>',
        TABIX   : '<li>Li, H., 2011. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27, 718-719. doi: <a href="https://doi.org/10.1093/bioinformatics/btq671">10.1093/bioinformatics/btq671</a></li>',
        GLIMPSE1: '<li>Rubinacci, S., Ribeiro, D.M., Hofmeister, R.J., Delaneau, O., 2021. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. Nat Genet 53, 120-126. doi: <a href="https://doi.org/10.1038/s41588-020-00756-0">10.1038/s41588-020-00756-0</a></li>',
        GLIMPSE2: '<li>Rubinacci, S., Hofmeister, R.J., Sousa da Mota, B., Delaneau, O., 2023. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. Nat Genet 55, 1088-1090. doi: <a href="https://doi.org/10.1038/s41588-023-01438-3">10.1038/s41588-023-01438-3</a></li>',
    ]

    if (steps.contains("all")) {
        steps = ["simulate", "panelprep", "impute", "validate"]
    }

    def reference_text = [
        tools.contains("beagle5")  ? tool_biblio.BEAGLE5  : "",
        steps.contains("panelprep") || steps.contains("validate") || steps.contains("simulate") || tools.contains("glimpse1") ? tool_biblio.SAM_BCFTOOLS : "",
        tools.contains("minimac4") ? tool_biblio.MINIMAC4 : "",
        tools.contains("stitch")   ? tool_biblio.STITCH   : "",
        tools.contains("quilt")    ? tool_biblio.QUILT    : "",
        tools.contains("quilt2")   ? tool_biblio.QUILT2   : "",
        tool_biblio.MULTIQC,
        steps.contains("panelprep") && compute_freq              ? tool_biblio.VCFLIB   : "",
        steps.contains("panelprep") && phase                     ? tool_biblio.SHAPEIT5 : "",
        steps.contains("validate") || tools.contains("glimpse1") ? tool_biblio.TABIX    : "",
        tools.contains("glimpse1") ? tool_biblio.GLIMPSE1 : "",
        tools.contains("glimpse2") ? tool_biblio.GLIMPSE2 : ""
    ].join(' ').trim().replaceAll("[,|.] +\\.", ".")

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml, steps, tools, normalize, remove_samples, compute_freq, phase) {
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
    meta["tool_citations"] = toolCitationText(
        steps, tools, normalize, remove_samples, compute_freq, phase
    ).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText(steps, tools, compute_freq, phase)

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
