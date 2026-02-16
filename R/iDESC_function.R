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
#' @param use_ind_zero Logical. If TRUE, adds a binary covariate to the zero-inflation component
#' to model the probability of structural zeros based on zero presence. Default is set to "TRUE"
#'
#' @return A list containing:
#' \itemize{
#'   \item{\strong{model_results}: A data frame containing the results of the differential expression analysis. For each gene, this includes the status of the model fit, estimated effect sizes (coefficients), dispersion parameters, zero-inflation estimates, statistical significance (raw p-values and FDR-corrected p-values), and model likelihood.}
#'   \item{\strong{total_cells}: A summary table containing the number of observations for each comparison group}
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
                covariates=NULL,cores=1, use_ind_zero=TRUE){
                
  if(!is.null(user_sf)){
    if(is.null(names(user_sf))){names(user_sf)<-colnames(mat)}else{user_sf<-user_sf[colnames(mat)]}
  }
    
  dat<-preprocessing(counts,meta,subject_var,group_var,sub_cell_filtering=5,gene_sub_filtering=0,gene_cell_filtering=0.05,ncell_filtering=1)
  mat<-dat$mat
  meta<-dat$meta
  norm_factor<-normalization_factors(mat,norm_opt,user_sf="SizeFactor")
  group<-meta[,group_var]
  subject<-meta[,subject_var]
  
  predict_pi<-zp_prediction(mat,norm_factor,span=0.05,loess_control="")
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
    "Alpha", group_coef_names, covar_coef_names,
    "Dispersion", "Sigma2", "Sigma2_ZI", "Theta",
    group_pval_names, covar_pval_names, "Minus2LogLik"
  )
  
  # Run model for each gene
  res_tb <- parallel::mclapply(1:nrow(filtered_mat),function(g){
  
    gene_name <- rownames(filtered_mat)[g]
    gene <- filtered_mat[g,]
    predict_pi_offset <- predict_pi[rownames(filtered_mat)[g]]
    tmp.df <- data.frame(y=gene,norm_sf=norm_factor,predict_pi_offset=predict_pi_offset,ind_zero=1*(gene==0),group=group,sub=as.numeric(factor(subject)))
  
    if (!is.null(covariates)) tmp.df[, covariates] <- meta[, covariates]
  
    # Create formulas for the model
    fixed_effects <- c("group", covariates)
    fixed_formula <- as.formula(paste("y ~", paste(fixed_effects, collapse = " + "), "+ offset(log(norm_sf)) + (1 | sub)"))
    if(use_ind_zero) {
      zi_formula <- ~ offset(predict_pi_offset) - 1 + ind_zero + (1 | sub)
    } else {
      zi_formula <- ~ offset(predict_pi_offset) - 1 + (1 | sub)
    }
  
    # Run model
    f1 <- tryCatch({
        suppressWarnings(suppressMessages(
            glmmTMB::glmmTMB(
              formula = fixed_formula,
              data = tmp.df,
              family = glmmTMB::nbinom2,
              ziformula = zi_formula
            )
          ))
      }, error = function(e) return(NULL))
  
    # Check convergence
    if (is.null(f1)) {
      fail_row <- data.frame(Gene = gene_name, Status = "model_crash", stringsAsFactors = FALSE)
      fail_row[expected_colnames] <- NaN
      return(return_fail("model_crash"))
    } else {
      status <- "OK"
    }
  
    summ <- summary(f1)
  
    est <-summ$coefficients$cond[, "Estimate"]
    pval <- summ$coefficients$cond[-1, "Pr(>|z|)"]
  
    if (is.null(pval)) {
      pval_names <- names(est)[-1]
      pval <- setNames(rep(NaN, length(pval_names)), pval_names)
    }
  
    zi_estimate <- NaN
    if (use_ind_zero) {
        zi_coefs <- summ$coefficients$zi
        if ("ind_zero" %in% rownames(zi_coefs)) {
            zi_estimate <- zi_coefs["ind_zero", "Estimate"]
        } else if (nrow(zi_coefs) > 0) {
            zi_estimate <- zi_coefs[1, "Estimate"]
        }
    }
  
    res1 <- c(est, 1 / summ$sigma, exp(f1$fit$par["theta"])^2, exp(f1$fit$par["thetazi"])^2,
              zi_estimate, pval, summ$AICtab["-2*log(L)"])
    
    names(res1) <- expected_colnames
    
    if (f1$fit$convergence != 0) status <- "convergence_error"
    if (any(is.na(summ$AICtab))) status <- "invalid_likelihood"
    if (any(is.na(summ$coefficients$cond[, "Pr(>|z|)"]))) status <- "invalid_pval"
  
    df <- data.frame(Gene = gene_name, Status = status, t(res1), stringsAsFactors = FALSE)
    return(df)
  },mc.cores=cores)
  
  res_tb <- plyr::rbind.fill(res_tb)
  
  # Include excluded genes (just expressed in one group) in the final table
  if (length(genes_to_exclude) > 0) {
    excluded_df <- data.frame(Gene = genes_to_exclude, Status = "one_group", stringsAsFactors = FALSE)
    excluded_df[, expected_colnames] <- NaN
    res_tb <- plyr::rbind.fill(res_tb, excluded_df)
  }
  
  rownames(res_tb) <- res_tb$Gene
  res_tb$Gene <- NULL
  
  # FDR correction
  pval_cols <- grep("^Pval_Beta_", colnames(res_tb), value = TRUE)
  for (pval_col in pval_cols) {
    fdr_col <- sub("^Pval_", "FDR_", pval_col)
    res_tb[[fdr_col]] <- NaN
    res_tb[[fdr_col]][which(res_tb$Status == "OK")] <- p.adjust(res_tb[[pval_col]][which(res_tb$Status == "OK")], method = "BH")
  }
  
  # Order FDR columns
  ordered_cols <- unlist(lapply(pval_cols, function(pval_col) {
    beta_col <- sub("^Pval_", "", pval_col)
    fdr_col <- sub("^Pval_", "FDR_", pval_col)
    c(beta_col, pval_col, fdr_col)
  }))
  
  final_col_order <- c(
    c("Status", "Alpha"),
    ordered_cols,
    setdiff(colnames(res_tb), c("Status", "Alpha", ordered_cols))
  )
  
  res_tb <- res_tb[, final_col_order]

  return(list(
    model_results = res_tb,
    total_cells = total_cells_per_group
  ))
}
