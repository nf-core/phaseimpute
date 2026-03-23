//
// Check if the contig names can be solved by adding/removing the `chr` prefix.
//
def diffChr(chr_target, chr_ref, file) {
    def diff = chr_ref - chr_target
    def prefix = (chr_ref - chr_target) =~ "chr" ? "chr" : "nochr"
    if (diff.size() != 0) {
        // Ensure that by adding/removing the prefix we can solve the problem
        def new_chr = []
        def to_rename = []
        if (prefix == "chr") {
            chr_target.each{ chr -> new_chr += "chr${chr}" }
            diff.each{ chr -> to_rename += chr.replace('chr', '') }
        } else {
            chr_target.each{ chr -> new_chr += chr.replace('chr', '') }
            diff.each{ chr -> to_rename += "chr${chr}" }
        }
        def new_diff = diff - new_chr
        if (new_diff.size() != 0) {
            def chr_names = new_diff.size() > params.max_chr_names ? new_diff[0..params.max_chr_names - 1] + ['...'] : new_diff
            error "Contig names: ${chr_names} absent from file: ${file} and cannot be solved by adding or removing the `chr` prefix."
        }
        diff = to_rename
    }
    return [diff, prefix]
}

//
// Check if the contig names in the input files match the reference contig names.
//
def checkChr(ch_chr, ch_input) {
    def chr_checked = ch_chr
        .combine(ch_input, by: 0)
        .map { meta, chr, file, index, lst ->
            [
                meta, file, index,
                chr.readLines()*.split(' ').collect { line -> line[0] },
                lst
            ]
        }
        .branch { meta, file, index, chr, lst ->
            def lst_diff = diffChr(chr, lst, file)
            def diff = lst_diff[0]
            def prefix = lst_diff[1]
            no_rename: diff.size() == 0
                return [meta, file, index]
            to_rename: true
                return [meta, file, index, diff, prefix]
        }
    return chr_checked
}
