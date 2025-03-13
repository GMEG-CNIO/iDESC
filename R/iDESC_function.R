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
#' @param cores Number of cores for parallel processing (default is 1)
#' @param loess_control Optional. If set to "control", the LOESS smoothing function will use
#' `loess.control(surface = "direct")`. This setting is recommended when working with sparse or
#' zero-heavy datasets, where the default LOESS method may struggle to define neighborhoods for
#' smoothing. Using "direct" computation ensures more stable and accurate smoothing but may be
#' slower. Activate this option if you encounter warnings like "k-d tree limited by memory" or
#' if the LOESS fit appears unstable.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\strong{model_results}: A data frame containing the results of the differential expression analysis. For each gene, this includes estimated effect sizes (coefficients), dispersion parameters, zero-inflation estimates, statistical significance (p-values), and model deviances.}
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
                cores=1){
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
  colnames(gene_group_expressed_df) <- levels(factor(group))

  # Total cells per group (only once)
  total_cells_per_group <- as.data.frame(as.list(table(group)))

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

  res_tb<-Reduce(plyr::rbind.fill,parallel::mclapply(1:nrow(filtered_mat),function(g){
    gene<-filtered_mat[g,]
    predict_pi_offset <- predict_pi[rownames(filtered_mat)[g]]
    tmp.df<-data.frame(y=gene,norm_sf=norm_factor,predict_pi_offset=predict_pi_offset,ind_zero=1*(gene==0),group=group,sub=as.numeric(factor(subject)))

    f1<-try(glmmTMB::glmmTMB(y ~ group + offset(log(norm_sf)) + (1 | sub), data = tmp.df, family = glmmTMB::nbinom2, zi = ~ offset(predict_pi_offset)-1+ind_zero+(1 | sub)))

    res1<-try(c(summary(f1)$coefficients$cond[,"Estimate"],1/summary(f1)$sigma,exp(f1$fit$par["theta"])^2,exp(f1$fit$par["thetazi"])^2,
                summary(f1)$coefficients$zi[1,"Estimate"],summary(f1)$coefficients$cond[-1,"Pr(>|z|)"],summary(f1)$AICtab["deviance"]))
    if('try-error' %in% class(res1)){
      res1<-c(Alpha=NA)
    }else{
      names(res1)<-c("Alpha",paste0("Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]), "Dispersion","Sigma2","Sigma2_ZI","Theta",
                     paste0("Pval_Beta_",names(summary(f1)$coefficients$cond[,"Estimate"])[-1]),"Deviance")

    }

    data.frame(t(res1))
  },mc.cores=cores))
  rownames(res_tb)<-rownames(filtered_mat)
  res_tb<-res_tb[which(!is.na(res_tb[,1])&!is.na(res_tb[,2])&!is.na(res_tb[,7])),]
  return(list(
    model_results = res_tb,
    problematic_genes_cells_expressed = problematic_genes_expressed,
    total_cells = total_cells_per_group
  ))
}
