#' iDESC: main function for DE analysis using zero-inflated negative binomial mixed model
#'
#' @param mat Count matrix
#' @param meta Data frame including information for cells
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta
#' @param norm_opt Option for normalizing factors
#' @param user_sf Option for user-specific normalizing factors
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject
#' @param gene_cell_filtering Filtering on genes based on expression across all cells
#' @param ncell_filtering Filtering on cells based on the number of genes expressed
#' @param span smoothing parameter for LOESS curve
#' @param loess_control Optional. If set to "control", the LOESS smoothing function will use
#' `loess.control(surface = "direct")`. This setting is recommended when working with sparse or
#' zero-heavy datasets, where the default LOESS method may struggle to define neighborhoods for
#' smoothing. Using "direct" computation ensures more stable and accurate smoothing but may be
#' slower. Activate this option if you encounter warnings like "k-d tree limited by memory" or
#' if the LOESS fit appears unstable.
#' @param covariates Optional. A vector of covariates to include in the model.
#' @param cores Number of cores for parallel processing (default is 1)
#'
#' @return A list containing:
#' \itemize{
#'   \item{\strong{model_results}: A data frame containing the results of the differential expression analysis. For each gene, this includes estimated effect sizes (coefficients), dispersion parameters, zero-inflation estimates, statistical significance (raw p-values and FDR-corrected p-values), and model deviances.}
#'   \item{\strong{problematic_genes}: A data frame listing genes that were excluded from the main analysis because they were expressed in cells from only one group. For each gene, the table reports the number of cells with non-zero expression in each group. This can be used for further diagnostic or exploratory analyses, such as assessing group-specific expression patterns.}
#' }
#' @export
#'
#' @examples
#'
#' library(iDESC)
#' data(IPF_example)
#' mat=IPF_example$mat
#' meta=IPF_example$meta
#' sequencing_depth=IPF_example$meta$sequencing_depth
#' result=iDESC(mat,meta,subject_var="subject",group_var="disease",norm_opt="User",user_sf = sequencing_depth,span = 0.7,cores = 4)

iDESC<-function(mat,meta,subject_var,group_var,norm_opt=c("SeqDepth","SizeFactor","User","None"),user_sf=NULL,
                sub_cell_filtering=5,gene_sub_filtering=0,gene_cell_filtering=0.05,ncell_filtering=1,span=0.05,loess_control="",
                covariates=NULL,cores=1){
  if(!is.null(user_sf)){
    if(is.null(names(user_sf))){names(user_sf)<-colnames(mat)}else{user_sf<-user_sf[colnames(mat)]}
  }
  dat<-preprocessing(mat,meta,subject_var,group_var,sub_cell_filtering=sub_cell_filtering,gene_sub_filtering=gene_sub_filtering,gene_cell_filtering=gene_cell_filtering,ncell_filtering=ncell_filtering)
  mat<-dat$mat
  meta<-dat$meta
  norm_factor<-normalization_factors(mat,norm_opt,user_sf)
  group<-meta[,group_var]
  subject<-meta[,subject_var]

  predict_pi<-zp_prediction(mat,norm_factor,span,loess_control)
  stopifnot("LOESS: zero-width neighborhood. make span bigger"=sum(is.na(predict_pi))==0)

  # Summary of expressed cells
  gene_group_expressed <- sapply(rownames(mat), function(gene_name) {
    gene_values <- mat[gene_name, ]
    tapply(gene_values > 0, group, sum)
  })
  gene_group_expressed_df <- as.data.frame(t(gene_group_expressed))
  gene_group_expressed_df <- gene_group_expressed_df[, colSums(!is.na(gene_group_expressed_df)) > 0]
  colnames(gene_group_expressed_df) <- levels(factor(group))

  # Total cells per group (only once)
  if (!is.factor(group)) {
    group <- factor(group)
  }
  group_clean <- droplevels(group)
  total_cells_per_group <- as.data.frame(as.list(table(group_clean)))
  rownames(total_cells_per_group) <- "Total_cells"

  # Identify and exclude genes expressed in only one group
  genes_to_exclude <- rownames(gene_group_expressed_df)[rowSums(gene_group_expressed_df > 0) < 2]

  # Filter summaries for problematic genes
  problematic_genes_expressed <- gene_group_expressed_df[genes_to_exclude, , drop = FALSE]

  # Issue warning if genes are excluded
  if (length(genes_to_exclude) > 0) {
    message("The following genes were excluded because they were expressed in only one group: ",
            paste(genes_to_exclude, collapse = ", "))
  }

  # Filter out problematic genes from the main matrix
  filtered_mat <- mat[!rownames(mat) %in% genes_to_exclude, ]

  # Prepare colnames for the results table when glmmTMB fails

  # Helper to get levels if factor, or just keep the variable name if numeric
  get_coef_names <- function(var_name, meta_df, prefix = "Beta_") {
    var_data <- meta_df[[var_name]]
    if (is.factor(var_data) || is.character(var_data)) {
      levels_var <- levels(factor(var_data))
      if (length(levels_var) <= 1) return(NULL)
      return(paste0(prefix, var_name, levels_var[-1]))  # drop reference level
    } else {
      return(paste0(prefix, var_name))
    }
  }

  # Get coefficients and p-values for group
  group_coef_names <- get_coef_names(group_var, meta, prefix = "Beta_")
  group_pval_names <- get_coef_names(group_var, meta, prefix = "Pval_Beta_")

  # For covariates
  covar_coef_names <- unlist(lapply(covariates, get_coef_names, meta_df = meta, prefix = "Beta_"))
  covar_pval_names <- unlist(lapply(covariates, get_coef_names, meta_df = meta, prefix = "Pval_Beta_"))

  # Combine all expected column names
  expected_colnames <- c(
    "Alpha",
    group_coef_names,
    covar_coef_names,
    "Dispersion", "Sigma2", "Sigma2_ZI", "Theta",
    group_pval_names,
    covar_pval_names,
    "Deviance"
  )

  convergence_failed_genes <- character(0)

  res_tb<-parallel::mclapply(1:nrow(filtered_mat),function(g){
    gene_name <- rownames(filtered_mat)[g]
    gene<-filtered_mat[g,]
    predict_pi_offset <- predict_pi[rownames(filtered_mat)[g]]
    tmp.df<-data.frame(y=gene,norm_sf=norm_factor,predict_pi_offset=predict_pi_offset,ind_zero=1*(gene==0),group=group,sub=as.numeric(factor(subject)))

    if (!is.null(covariates)) {
      tmp.df[, covariates] <- meta[, covariates]
    }

    # Create formulas for the model
    fixed_effects <- c("group", covariates)
    fixed_formula <- as.formula(paste("y ~", paste(fixed_effects, collapse = " + "), "+ offset(log(norm_sf)) + (1 | sub)"))
    zi_formula <- ~ offset(predict_pi_offset) - 1 + ind_zero + (1 | sub)

    f1 <- tryCatch(
      {
        suppressWarnings(
          suppressMessages(
            glmmTMB::glmmTMB(
              formula = fixed_formula,
              data = tmp.df,
              family = glmmTMB::nbinom2,
              ziformula = zi_formula
            )
          )
        )
      },
      error = function(e) {
        convergence_failed_genes <- c(convergence_failed_genes, gene_name)
        message(paste0("Model convergence failed for gene: ", gene_name))
        return(NULL)
      }
    )

    # Check convergence
    if (is.null(f1)) return(NULL)

    est <- summary(f1)$coefficients$cond[, "Estimate"]
    pval <- summary(f1)$coefficients$cond[-1, "Pr(>|z|)"]
    if (is.null(pval)) {
      pval_names <- names(est)[-1]
      pval <- setNames(rep(NA, length(pval_names)), pval_names)
    }
    res1 <- c(est, 1 / summary(f1)$sigma, exp(f1$fit$par["theta"])^2, exp(f1$fit$par["thetazi"])^2,
              summary(f1)$coefficients$zi[1, "Estimate"], pval, summary(f1)$AICtab["deviance"])
    names(res1) <- c("Alpha",
                      paste0("Beta_", names(est)[-1]),
                      "Dispersion", "Sigma2", "Sigma2_ZI", "Theta",
                      paste0("Pval_Beta_", names(est)[-1]),
                      "Deviance")

    df <- data.frame(Gene = gene_name, t(res1), stringsAsFactors = FALSE)
    return(df)
  },mc.cores=cores)

  res_tb <- plyr::rbind.fill(res_tb)
  rownames(res_tb) <- res_tb$Gene
  res_tb$Gene <- NULL

  # Add convergence-failure genes to the problematic table
  if (length(convergence_failed_genes) > 0) {
    convergence_failed_expr <- gene_group_expressed_df[convergence_failed_genes, , drop = FALSE]

    convergence_failed_expr$Reason <- "model_convergence"
    convergence_failed_expr$Gene <- rownames(convergence_failed_expr)

    convergence_failed_df <- convergence_failed_expr
  } else {
    convergence_failed_df <- NULL
  }

  # Combine both problematic gene sets into one
  if (nrow(problematic_genes_expressed) > 0) {
    problematic_expr <- problematic_genes_expressed
    problematic_expr$Reason <- "one_group"
    problematic_expr$Gene <- rownames(problematic_expr)
  } else {
    problematic_expr <- NULL
  }

  problematic_combined <- plyr::rbind.fill(problematic_expr, convergence_failed_df)
  if (!is.null(problematic_combined)) {
    problematic_combined <- problematic_combined[order(problematic_combined$Gene), ]
    problematic_combined <- problematic_combined[, !(colnames(problematic_combined) %in% "Gene")]
  }

  # FDR correction
  pval_cols <- grep("^Pval_Beta_", colnames(res_tb), value = TRUE)
  for (pval_col in pval_cols) {
    fdr_col <- sub("^Pval_", "FDR_", pval_col)
    res_tb[[fdr_col]] <- p.adjust(res_tb[[pval_col]], method = "BH")
  }

  # Order FDR columns
  ordered_cols <- unlist(lapply(pval_cols, function(pval_col) {
    fdr_col <- sub("^Pval_", "FDR_", pval_col)
    c(pval_col, fdr_col)
  }))
  final_col_order <- c(
    setdiff(colnames(res_tb), c(ordered_cols, "Deviance")),
    ordered_cols,
    "Deviance"
  )

  return(list(
    model_results = res_tb,
    problematic_genes_cells_expressed = problematic_combined,
    total_cells = total_cells_per_group
  ))
}
