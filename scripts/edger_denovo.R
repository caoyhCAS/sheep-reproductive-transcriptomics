#!/usr/bin/env Rscript
# Method-based reconstruction for DOI 10.1038/s41437-018-0090-1.
# These are not the authors' original scripts. edgeR's study version is unknown.

abort <- function(message) stop(message, call. = FALSE)

parse_options <- function(args, value_names, required_names) {
  flags <- c("allow-pseudoreplication", "validate-only", "help")
  result <- list()
  i <- 1L
  while (i <= length(args)) {
    name <- sub("^--", "", args[i])
    if (!grepl("^--", args[i]) || !name %in% c(value_names, flags))
      abort(paste("Unknown option:", args[i]))
    if (!is.null(result[[name]])) abort(paste("Duplicate option:", name))
    if (name %in% flags) {
      result[[name]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args) || grepl("^--", args[i + 1L]))
        abort(paste("Missing value for:", name))
      value <- args[i + 1L]
      if (!nzchar(trimws(value)) || grepl("[\r\n\t]", value))
        abort(paste("Invalid empty/multiline option:", name))
      result[[name]] <- value
      i <- i + 2L
    }
  }
  if (isTRUE(result$help)) return(result)
  missing <- required_names[!required_names %in% names(result)]
  if (length(missing)) abort(paste("Missing required options:", paste(missing, collapse = ", ")))
  if (identical(result$numerator, result$denominator)) abort("numerator and denominator must differ")
  result
}

read_tsv <- function(path) {
  if (!file.exists(path) || dir.exists(path)) abort(paste("Input file not found:", path))
  table <- utils::read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                            check.names = FALSE, stringsAsFactors = FALSE,
                            colClasses = "character", na.strings = character())
  if (anyDuplicated(names(table))) abort(paste("Duplicate column names:", path))
  if (!nrow(table)) abort(paste("Empty input table:", path))
  table
}

check_ids <- function(ids, label) {
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
    abort(paste(label, "must be nonempty and unique"))
  if (any(ids != trimws(ids)) || any(grepl("[\r\n\t]", ids)))
    abort(paste(label, "must not contain padding or tab/newline characters"))
}

read_comparison <- function(options, integer_counts = FALSE) {
  input <- read_tsv(options$counts)
  if (ncol(input) < 3L || names(input)[1L] != "gene_id")
    abort("Counts require first column gene_id and at least two sample columns")
  check_ids(input$gene_id, "gene_id values")
  check_ids(names(input)[-1L], "Count sample names")
  samples <- read_tsv(options$samples)
  if (!all(c("sample_id", "condition", "individual_id") %in% names(samples)))
    abort("Samples require sample_id, condition, individual_id columns")
  check_ids(samples$sample_id, "Metadata sample_id values")
  if (any(!nzchar(trimws(samples$condition))) || any(!nzchar(trimws(samples$individual_id))))
    abort("condition and individual_id must not be empty")
  for (field in c("condition", "individual_id")) {
    if (any(samples[[field]] != trimws(samples[[field]])) || any(grepl("[\r\n\t]", samples[[field]])))
      abort(paste(field, "must not contain padding or tab/newline characters"))
  }
  if (!setequal(names(input)[-1L], samples$sample_id))
    abort("Count sample columns and metadata sample_id values must match exactly")
  samples <- samples[match(names(input)[-1L], samples$sample_id), , drop = FALSE]
  values <- suppressWarnings(as.numeric(as.matrix(input[, -1L, drop = FALSE])))
  counts <- matrix(values, nrow = nrow(input), ncol = ncol(input) - 1L,
                   dimnames = list(input$gene_id, samples$sample_id))
  if (any(!is.finite(counts)) || any(counts < 0))
    abort("Counts must be finite, nonnegative numbers; no missing values or expression units")
  if (integer_counts && (any(counts != floor(counts)) || any(counts > .Machine$integer.max)))
    abort("DESeq requires integer read counts within R's integer range; no rounding, TPM, FPKM or fractional counts")
  selected <- samples$condition %in% c(options$denominator, options$numerator)
  samples <- samples[selected, , drop = FALSE]
  counts <- counts[, selected, drop = FALSE]
  if (!all(c(options$numerator, options$denominator) %in% samples$condition))
    abort("Both numerator and denominator must have samples")
  if (any(!is.finite(colSums(counts))) || any(colSums(counts) <= 0))
    abort("Every selected library must have a positive finite count total")
  groups <- factor(samples$condition, levels = c(options$denominator, options$numerator))
  biological_n <- vapply(split(samples$individual_id, groups), function(x) length(unique(x)), integer(1))
  repeated <- anyDuplicated(samples$individual_id) > 0L
  nonindependent <- repeated || any(biological_n < 2L)
  if (nonindependent && !isTRUE(options[["allow-pseudoreplication"]]))
    abort(paste("Repeated individuals or fewer than two biological individuals in a group:",
                "--allow-pseudoreplication is required for a descriptive historical comparison"))
  if (nonindependent)
    message("DESCRIPTIVE ONLY: repeated libraries/individuals do not supply independent biological replication; nominal P values cannot support population inference.")
  list(counts = counts, samples = samples, groups = groups, biological_n = biological_n,
       nonindependent = nonindependent)
}

output_paths <- function(output) c(output, paste0(output, ".provenance.tsv"), paste0(output, ".sessionInfo.txt"))

check_output_paths <- function(output) {
  paths <- output_paths(output)
  exists <- vapply(paths, function(path) file.exists(path) || isTRUE(nzchar(Sys.readlink(path))), logical(1))
  if (any(exists)) abort(paste("Refusing to overwrite:", paste(paths[exists], collapse = ", ")))
}

new_results <- function(data, options) {
  data.frame(gene_id = rownames(data$counts), log2fc = NA_real_, padj = NA_real_,
             status = "NOTEST", numerator = options$numerator, denominator = options$denominator,
             pvalue = NA_real_, mean_expression = NA_real_,
             analysis_scope = if (data$nonindependent) "descriptive_nonindependent_libraries" else "unpaired_two_group",
             stringsAsFactors = FALSE)
}

mark_result_status <- function(result) {
  valid <- !is.na(result$log2fc) & is.finite(result$padj) & result$padj >= 0 & result$padj <= 1
  # +/-Inf log2 fold changes can represent expression in only one condition.
  result$status <- ifelse(valid, "OK", "NOTEST")
  result$padj[!valid] <- NA_real_
  result
}

write_results <- function(result, data, options, package, version, details) {
  check_output_paths(options$output)
  parent <- dirname(options$output)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) abort("Unable to create output directory")
  details <- c(
    method_status = "Method-based reconstruction; not original author code",
    package = package, installed_version = version,
    numerator = options$numerator, denominator = options$denominator,
    log2fc_orientation = "numerator / denominator",
    selected_samples = paste(data$samples$sample_id, collapse = ","),
    selected_individuals = paste(data$samples$individual_id, collapse = ","),
    biological_individuals_per_group = paste(names(data$biological_n), data$biological_n, sep = ":", collapse = ","),
    pseudoreplication_override = as.character(isTRUE(options[["allow-pseudoreplication"]])),
    descriptive_only = as.character(data$nonindependent),
    counts_md5 = unname(tools::md5sum(options$counts)),
    samples_md5 = unname(tools::md5sum(options$samples)),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), details)
  utils::write.table(result, file = options$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  utils::write.table(data.frame(key = names(details), value = unname(details)),
                     file = paste0(options$output, ".provenance.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), paste0(options$output, ".sessionInfo.txt"))
}

run_edger <- function(args = commandArgs(trailingOnly = TRUE)) {
  required <- c("counts", "samples", "numerator", "denominator", "output")
  options <- parse_options(args, c(required, "dispersion", "expected-version"), required)
  if (isTRUE(options$help)) {
    cat("Method reconstruction, not original code.\n",
        "Rscript scripts/edger_denovo.R --counts counts.tsv --samples contrast.tsv\n",
        "  --numerator ovary --denominator endometrium --output result.tsv\n",
        "  [--allow-pseudoreplication] [--dispersion NUMBER] [--expected-version VERSION] [--validate-only]\n", sep = "")
    return(invisible(NULL))
  }
  # Validation-only deliberately works with base R without statistical packages.
  if (!isTRUE(options[["validate-only"]]) && !requireNamespace("edgeR", quietly = TRUE))
    abort("edgeR is not installed; no package is installed or substituted automatically")
  data <- read_comparison(options, integer_counts = FALSE)
  check_output_paths(options$output)
  dispersion <- NULL
  if (!is.null(options$dispersion)) {
    dispersion <- suppressWarnings(as.numeric(options$dispersion))
    if (length(dispersion) != 1L || !is.finite(dispersion) || dispersion <= 0)
      abort("dispersion must be a positive finite caller-supplied value")
  }
  if (ncol(data$counts) <= 2L && is.null(dispersion))
    abort("No residual library replication: provide an explicitly justified --dispersion; no biological variance is invented")
  tested <- rowSums(data$counts) > 0
  if (!any(tested)) abort("No nonzero features available")
  if (isTRUE(options[["validate-only"]])) {
    cat("VALIDATION ONLY: input contracts accepted; edgeR availability, version and numerical results not tested.\n")
    return(invisible(NULL))
  }
  version <- as.character(utils::packageVersion("edgeR"))
  if (!is.null(options[["expected-version"]]) && version != options[["expected-version"]])
    abort(paste("edgeR version mismatch:", version, "!=", options[["expected-version"]]))
  message("edgeR ", version, "; study version unreported; classical TMM/exactTest choices are reconstruction assumptions.")
  y <- edgeR::DGEList(counts = data$counts[tested, , drop = FALSE], group = data$groups)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  if (is.null(dispersion)) {
    y <- edgeR::estimateCommonDisp(y)
    y <- edgeR::estimateTagwiseDisp(y)
    selected_dispersion <- "tagwise"
  } else selected_dispersion <- dispersion
  fit <- edgeR::exactTest(y, pair = c(options$denominator, options$numerator),
                         dispersion = selected_dispersion, rejection.region = "doubletail",
                         big.count = 900, prior.count = 0.125)
  table <- fit$table
  result <- new_results(data, options)
  index <- match(rownames(table), result$gene_id)
  result$log2fc[index] <- table$logFC
  result$pvalue[index] <- table$PValue
  result$padj[index] <- stats::p.adjust(table$PValue, method = "BH")
  result$mean_expression[index] <- rowMeans(edgeR::cpm(y, normalized.lib.sizes = TRUE))
  result <- mark_result_status(result)
  write_results(result, data, options, "edgeR", version, c(
    reported_version = "unreported", normalization = "TMM", test = "classical exactTest; doubletail; prior.count=0.125; big.count=900",
    dispersion = if (is.null(dispersion)) "estimateCommonDisp + estimateTagwiseDisp; installed defaults" else as.character(dispersion),
    feature_filter = "all-zero rows untested; no additional low-count filter", mean_expression_units = "normalized counts per million",
    paper_threshold = "abs(log2fc)>=2 and padj<=0.01; apply downstream",
    design_limit = "individual_id audits independence only; it is not included as a blocking factor"))
  invisible(result)
}

if (sys.nframe() == 0L) {
  tryCatch(run_edger(), error = function(e) { message("ERROR: ", conditionMessage(e)); quit(status = 2L) })
}
