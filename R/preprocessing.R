#' preprocessing: filtering step on genes/cells
#'
#' @param mat  Count matrix
#' @param meta Data frame including information for cells
#' @param subject_var The name of subject information in meta
#' @param group_var The name of group/disease information for DE analysis in meta
#' @param sub_cell_filtering Filtering on cells within each subject
#' @param gene_sub_filtering Filtering on genes based on expression in each subject
#' @param gene_cell_filtering Filtering on genes based on expression across all cells
#' @param ncell_filtering Filtering on cells based on the number of genes expressed
#' @param cores Number of cores for parallel processing (default is 1)
#'
#' @return A list of processed count matrix and data frame
#' @export

preprocessing<-function(mat,meta,subject_var,group_var,sub_cell_filtering,gene_sub_filtering,gene_cell_filtering,ncell_filtering,cores=1){
  #meta is a df with rownames=cell_barcode, columns includes:subject,group
  #sub_cell_filtering is an integer
  stopifnot("Gene filtering should be based on an indicated proportion"=(gene_sub_filtering>=0&gene_sub_filtering<=1))
  stopifnot("Gene filtering should be based on an indicated proportion"=(gene_cell_filtering>=0&gene_cell_filtering<=1))
  stopifnot("Cell barcodes in meta data and matrix should be the same"=(colnames(mat)==rownames(meta)))

  mat<-as.matrix(mat)

  #Remove subject with <"sub_cell_filtering" cells
  barcode_orig<-rownames(meta)
  subject_name<-names(table(meta[,subject_var]))
  remain_sub<-subject_name[which(table(meta[,subject_var])>=sub_cell_filtering)]
  sub_cell<-barcode_orig[which(meta[,subject_var]%in%remain_sub)]
  expressed_cell<-barcode_orig[which(Matrix::colSums(mat!=0)>0)]
  remain_cell<-intersect(sub_cell,expressed_cell)
  meta<-meta[remain_cell,]

  #Keep genes expressed at least a% cells in b% subject in either group
  if(length(levels(factor(meta[,group_var])))>3){
    print('"Group" has more than 3 levels and was treated as continuous.')
    group_sub_list<-list(subject_name)
  }else{
    group_sub_list<-parallel::mclapply(1:length(levels(factor(meta[,group_var]))),function(i){
      unique(meta[which(meta[,group_var]==levels(factor(meta[,group_var]))[i]),subject_var])
    },mc.cores=cores)
  }

  #sub_prop is a gene*sub matrix
  sub_prop_list<-parallel::mclapply(remain_sub,function(sub){
    rowMeans(mat[,barcode_orig[which(meta[,subject_var]==sub)]]!=0,na.rm=T)
  },mc.cores=cores)

  # list to matrix
  sub_prop <- do.call(cbind, sub_prop_list)
  colnames(sub_prop) <- remain_sub

  remain_gene<-rownames(mat)[which(rowSums(sapply(group_sub_list,function(sub_list){
    rowSums(sub_prop[,as.character(sub_list)]>=gene_cell_filtering,na.rm=T)>=(gene_sub_filtering*length(sub_list))
  }))>0)]

  # filter cells by number of expressed genes (ncell_filtering)
  mat_filtered <- mat[remain_gene, remain_cell]
  n_genes_per_cell <- Matrix::colSums(mat_filtered != 0)
  final_cells <- names(n_genes_per_cell[n_genes_per_cell >= ncell_filtering])

  # Subset matrix and metadata to keep only those cells
  mat_filtered <- mat_filtered[, final_cells]
  meta_filtered <- meta[final_cells, ]

  return(list(mat = mat_filtered, meta = meta_filtered))
}
