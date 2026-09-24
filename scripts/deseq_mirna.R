#!/usr/bin/env Rscript
# Method reconstruction for DOI 10.1038/s41437-018-0090-1; DESeq v1.22.0 only.
# Source shared input/output contracts without running the edgeR entry point.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript so the script directory can be resolved", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
source(file.path(script_dir, "edger_denovo.R"), local = TRUE)

run_deseq <- function(args = commandArgs(trailingOnly = TRUE)) {
  required <- c("counts", "samples", "features", "numerator", "denominator", "output",
                "dispersion-method", "sharing-mode", "fit-type")
  options <- parse_options(args, required, required)
  if (isTRUE(options$help)) {
    cat("Method reconstruction, not original code. DESeq v1.22.0 only; DESeq2 is not substituted.\n",
        "Rscript scripts/deseq_mirna.R --counts counts.tsv --samples contrast.tsv --features features.tsv\n",
        "  --numerator mouflon --denominator finnsheep --output result.tsv\n",
        "  --dispersion-method pooled|per-condition|blind --sharing-mode maximum|fit-only|gene-est-only\n",
        "  --fit-type parametric|local [--allow-pseudoreplication] [--validate-only]\n", sep = "")
    return(invisible(NULL))
  }
  if (!options[["dispersion-method"]] %in% c("pooled", "per-condition", "blind"))
    abort("dispersion-method must be pooled, per-condition or blind")
  if (!options[["sharing-mode"]] %in% c("maximum", "fit-only", "gene-est-only"))
    abort("sharing-mode must be maximum, fit-only or gene-est-only")
  if (!options[["fit-type"]] %in% c("parametric", "local")) abort("fit-type must be parametric or local")
  # Fail before processing data when analysis was requested with an absent/wrong package.
  if (!isTRUE(options[["validate-only"]])) {
    if (!requireNamespace("DESeq", quietly = TRUE)) abort("DESeq v1.22.0 is not installed; DESeq2 is not a substitute")
    version <- as.character(utils::packageVersion("DESeq"))
    if (version != "1.22.0") abort(paste("DESeq v1.22.0 required; installed", version, "is not silently substituted"))
  }
  data <- read_comparison(options, integer_counts = TRUE)
  check_output_paths(options$output)
  features <- read_tsv(options$features)
  if (!all(c("gene_id", "class") %in% names(features))) abort("Features require gene_id and class columns")
  check_ids(features$gene_id, "Feature gene_id values")
  if (!setequal(features$gene_id, rownames(data$counts))) abort("Feature IDs must exactly match count gene_id values")
  features <- features[match(rownames(data$counts), features$gene_id), , drop = FALSE]
  if (any(!features$class %in% c("sheep", "conserved", "novel")))
    abort("miRNA class must be exactly sheep, conserved or novel")
  tested <- features$class %in% c("sheep", "conserved") & rowSums(data$counts) > 0
  if (!any(tested)) abort("No nonzero sheep/conserved miRNAs available for testing")
  if (any(colSums(data$counts[tested, , drop = FALSE]) <= 0)) abort("Each library requires positive sheep/conserved counts")
  if (!any(rowSums(data$counts[tested, , drop = FALSE] > 0) == ncol(data$counts)))
    abort("Legacy size-factor estimation needs at least one retained feature positive in every library; no pseudocount is added")
  if (all(table(data$groups) < 2L) && options[["dispersion-method"]] != "blind")
    abort("No within-condition library replication: only explicitly selected blind dispersion is supported")
  if (isTRUE(options[["validate-only"]])) {
    cat("VALIDATION ONLY: input contracts accepted; DESeq availability, version and numerical results not tested.\n")
    return(invisible(NULL))
  }
  count_data <- data$counts[tested, , drop = FALSE]
  storage.mode(count_data) <- "integer"
  cds <- DESeq::newCountDataSet(count_data, data$groups)
  cds <- DESeq::estimateSizeFactors(cds)
  cds <- DESeq::estimateDispersions(cds, method = options[["dispersion-method"]],
                                   sharingMode = options[["sharing-mode"]], fitType = options[["fit-type"]])
  # DESeq reports condB / condA; always use the requested denominator first.
  table <- DESeq::nbinomTest(cds, condA = options$denominator, condB = options$numerator)
  result <- new_results(data, options)
  result$class <- features$class
  index <- match(table$id, result$gene_id)
  result$log2fc[index] <- table$log2FoldChange
  result$pvalue[index] <- table$pval
  result$padj[index] <- table$padj
  result$mean_expression[index] <- table$baseMean
  result <- mark_result_status(result)
  write_results(result, data, options, "DESeq", version, c(
    reported_version = "1.22.0", normalization = "DESeq estimateSizeFactors; no transformed expression input",
    dispersion_method = options[["dispersion-method"]], sharing_mode = options[["sharing-mode"]], fit_type = options[["fit-type"]],
    test = "nbinomTest; Benjamini-Hochberg adjusted P values",
    feature_filter = "only nonzero sheep/conserved rows modeled; novel/all-zero rows are NOTEST",
    features_md5 = unname(tools::md5sum(options$features)), mean_expression_units = "DESeq normalized base mean",
    paper_threshold = "abs(log2fc)>=1 and padj<=0.05; apply downstream",
    reconstruction_assumptions = "dispersion/fit/sharing choices and normalization feature universe were not reported",
    design_limit = "individual_id audits independence only; it is not included as a blocking factor"))
  invisible(result)
}

if (sys.nframe() == 0L) {
  tryCatch(run_deseq(), error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 2L) })
}
