suppressPackageStartupMessages({library(Seurat);library(Matrix);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
out <- normalizePath(args[1]);revision <- file.path(out,'figure_revision_20260922')
counts <- as(readMM(file.path(revision,'marker_counts.mtx')),'CsparseMatrix')
genes <- read.csv(file.path(revision,'marker_genes.csv'))$gene
meta <- read.csv(file.path(revision,'marker_metadata_PRIVATE.csv'),check.names=FALSE)
stopifnot(ncol(counts)==42806,nrow(counts)==length(genes),all(c('GFP','mCherry') %in% genes),all(meta$total_counts>0))
dimnames(counts) <- list(genes,meta$cell)
obj <- CreateSeuratObject(counts=counts,min.cells=0,min.features=0)
# Keep original full-library size factors; never normalize the marker subset itself.
normalized <- counts %*% Diagonal(x=10000/meta$total_counts)
normalized@x <- log1p(normalized@x)
dimnames(normalized) <- dimnames(counts)
obj <- SetAssayData(obj,assay='RNA',layer='data',new.data=normalized)
states <- read.csv(file.path(out,'cluster_state_inference','state_summary.csv'))
stopifnot(identical(states$rank,1:19),sum(states$n_cells)==42806)
display <- setNames(paste0(states$state,' (n=',format(states$n_cells,big.mark=',',trim=TRUE),')'),sub('cl','',states$original_cluster))
Idents(obj) <- factor(unname(display[as.character(meta$cluster)]),levels=rev(unname(display)))
p <- DotPlot(obj,features=genes,cols=c('lightgrey','#08306b'),dot.scale=7,
             scale=TRUE,cluster.idents=FALSE,col.min=-2.5,col.max=2.5,scale.by='radius') +
     RotatedAxis() + theme_bw(base_size=14) +
     theme(axis.text.x=element_text(size=14,angle=90,hjust=1,vjust=.5,face='italic'),
           axis.text.y=element_text(size=13),axis.title=element_blank(),
           panel.grid.major=element_line(color='grey92',linewidth=.3),
           panel.grid.minor=element_blank(),legend.position='right')
write.csv(p$data,file.path(revision,'seurat_dotplot_data.csv'),row.names=FALSE)
ggsave(file.path(out,'candidate_r1.1_markers.png'),p,width=19,height=10,dpi=180)
ggsave(file.path(revision,'candidate_r1.1_markers_seurat.pdf'),p,width=19,height=10)
cat('Seurat version:',as.character(packageVersion('Seurat')),'\n')
cat('PASS: actual Seurat DotPlot; 19 ordered putative states,',length(genes),'genes, including GFP and mCherry.\n')
