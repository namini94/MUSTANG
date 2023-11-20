# Summary

This tutorial demonstrates an analysis workflow of MUSTANG on a mouse brain multi-sample 10X Visium ST data.
It includes loading the dataset, pre-processing, concatenating, spot similarity graph construction and visualization. 

# Data

Here in this tutorial, we focus on the multi-section 10X Visium mouse brain spatial transcriptomics dataset:

- [Mouse brain section 1 (Sagittal-Anterior)
dataset](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-0-0)
- [Mouse brain section 1 (Sagittal-Posterior)
dataset](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-0-0)
- [Mouse brain section 2 (Sagittal-Anterior)
dataset](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard-1-0-0)
- [Mouse brain section 2 (Sagittal-Posterior)
dataset](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-posterior-1-standard-1-0-0)


In particular, we are interested in the filtered count matrix and the
spatial positions of the barcodes. We can download these files into a
folder, which we will use to load them into a
[`SpatialExperiment`](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html)
object.

# Loading Data

Here, in this section we load and combine ST sections to construct the spots transcriptional graph.

``` r
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(harmony)
  library(BayesSpace)
})
set.seed(100)

matrix_dir="~/Documents/Spatial/mouse_brain/Sec1_Anterior/outs/"

sec1_anterior <- SpatialExperiment::read10xVisium(samples = matrix_dir,
                                                  type = "HDF5",
                                                  data = "filtered",load = T)

colData(sec1_anterior)$sample_id<-"Sec1_Anterior"

matrix_dir="~/Documents/Spatial/mouse_brain/Sec1_Posterior/outs/"

sec1_posterior <- SpatialExperiment::read10xVisium(samples = matrix_dir,
                                                  type = "HDF5",
                                                  data = "filtered",load = T)

colData(sec1_posterior)$sample_id<-"Sec1_Posterior"

matrix_dir="~/Documents/Spatial/mouse_brain/Sec2_Anterior/outs/"

sec2_anterior <- SpatialExperiment::read10xVisium(samples = matrix_dir,
                                                  type = "HDF5",
                                                  data = "filtered",load = T)

colData(sec2_anterior)$sample_id<-"Sec2_Anterior"

matrix_dir="~/Documents/Spatial/mouse_brain/Sec2_Posterior/outs/"

sec2_posterior <- SpatialExperiment::read10xVisium(samples = matrix_dir,
                                                   type = "HDF5",
                                                   data = "filtered",load = T)

colData(sec2_posterior)$sample_id<-"Sec2_Posterior"

rowData(sec1_anterior)$is.HVG = NULL 
rowData(sec1_posterior)$is.HVG = NULL 
rowData(sec2_anterior)$is.HVG = NULL 
rowData(sec2_posterior)$is.HVG = NULL 

for(i in 1:nrow(colData(sec1_anterior))){
  colData(sec1_anterior)@rownames[i]<-paste0("Sec1_Ant_",colData(sec1_anterior)@rownames[i])
}
for(i in 1:nrow(colData(sec1_posterior))){
  colData(sec1_posterior)@rownames[i]<-paste0("Sec1_Post_",colData(sec1_posterior)@rownames[i])
}
for(i in 1:nrow(colData(sec2_anterior))){
  colData(sec2_anterior)@rownames[i]<-paste0("Sec2_Ant_",colData(sec2_anterior)@rownames[i])
}
for(i in 1:nrow(colData(sec2_posterior))){
  colData(sec2_posterior)@rownames[i]<-paste0("Sec2_Post_",colData(sec2_posterior)@rownames[i])
}

#Combine into 1 SCE and preprocess
sce.combined = cbind(sec2_anterior, sec1_anterior, sec2_posterior, sec1_posterior, deparse.level = 1)
sce.combined = spatialPreprocess(sce.combined, n.PCs = 50, n.HVGs=2000,assay.type="logcounts") #lognormalize, PCA

sce.combined = runUMAP(sce.combined, dimred = "PCA")
colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")

sce.combined

```


    ## class: SpatialExperiment 
    ## dim: 31053 12167  
    ## metadata(1): BayesSpace.data
    ## assays(2): counts logcounts
    ## rownames(31053): ENSMUSG00000051951 ENSMUSG00000089699 ... ENSMUSG00000096730
    ##   ENSMUSG00000095742
    ## rowData names(2): symbol is.HVG
    ## colnames(12167): Sec2_Ant_AAACAAGTATCTCCCA-1 Sec2_Ant_AAACACCAATAACTGC-1 ...
    ## Sec1_Post_TTGTTTCATTAGTCTA-1 Sec1_Post_TTGTTTCCATACAACT-1
    ## colData names(5): in_tissue array_row array_col sample_id sizeFactor
    ## reducedDimNames(2): PCA UMAP
    ## mainExpName: NULL
    ## altExpNames(0):
    ## spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
    ## imgData names(4): sample_id image_id data scaleFactor

   
# Batch Correction
Correcting possible batch effects when performing multi-sample analysis is necessary. Here, we perform batch correction and visualize the UMAP embedding of spots before and after batch correction:


``` r
sce.combined = runUMAP(sce.combined, dimred = "PCA")
colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP"))) +
  geom_point(alpha=0.3,aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_id)),size=0.7) +
  labs(color = "Sample") +
  theme_classic() +
  scale_color_manual(breaks = c("Sec1_Anterior", "Sec1_Posterior", "Sec2_Anterior","Sec2_Posterior"),
                     values=c("steelblue1", "peru", "green","mediumpurple"))+
  theme(axis.line = element_line(colour = 'black', size = 1.5))+
  coord_fixed()

```
<img src="https://github.com/namini94/MUSTANG/blob/main/Miscel/Mouse_Brain_Markdown_Figs/MB_NoBatch.png" width="50%" height="50%">
Now, we do batch correction and visualize the umap embedding again:

``` r
colData(sce.combined)$sample_id<-as.factor(colData(sce.combined)$sample_id)

sce.combined = RunHarmony(sce.combined, c("sample_id"), verbose = T)
sce.combined = runUMAP(sce.combined, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce.combined, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY"))) +
  geom_point(alpha=0.3,aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_id)),size=0.7) +
  labs(color = "Sample") +
  theme_classic()+
  scale_color_manual(breaks = c("Sec1_Anterior", "Sec1_Posterior", "Sec2_Anterior","Sec2_Posterior"),
                     values=c("steelblue1", "peru", "green","mediumpurple"))+
  theme(axis.line = element_line(colour = 'black', size = 1.5))+
  coord_fixed()
```
<img src="https://github.com/namini94/MUSTANG/blob/main/Miscel/Mouse_Brain_Markdown_Figs/MB_wBatchCorrection.png" width="50%" height="50%">

# KNN Graph & Louvain Clustering

``` r
harmony<-data.frame(reducedDim(sce.combined, "HARMONY"))
harmony_umap<-data.frame(reducedDim(sce.combined, "UMAP.HARMONY"))
k <- 50
tempcom <- MERINGUE::getClusters(harmony, k, weight=TRUE, method = igraph::cluster_louvain)

dat <- data.frame("emb1" = harmony_umap$UMAP1,
                  "emb2" = harmony_umap$UMAP2,
                  "Cluster" = tempCom)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(alpha=0.4,ggplot2::aes(x = emb1, y = emb2,
                                             color = Cluster), size = 0.9) +
  
  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +
  
  ggplot2::labs(title = "",
                x = "UMAP1",
                y = "UMAP2") +
  
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black"),
                 axis.text.y = ggplot2::element_text(color = "black"),
                 axis.title.y = ggplot2::element_text(),
                 axis.title.x = ggplot2::element_text(),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text( colour = "black"),
                 legend.title = ggplot2::element_text( colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1.5, colour = "black")
                 # legend.position="none"
  ) +
  
  
  ggplot2::coord_fixed()

plt

```
<img src="https://github.com/namini94/MUSTANG/blob/main/Miscel/Mouse_Brain_Markdown_Figs/Louvain_Clusters.png" width="50%" height="50%">

# Spots Transcriptional Graph
Now that we have identified the spots trancriptional clusters, we can construct and store the edges of spots transcriptional graph.

``` r
sample_ID_Sec2Ant<-matrix(0,nrow(colData(sec2_anterior)),1)
for(i in 1:nrow(colData(sec2_anterior))){
  sample_ID_Sec2Ant[i,1]<-c("Sec2_Anterior")
}
sample_ID_Sec1Ant<-matrix(0,nrow(colData(sec1_anterior)),1)
for(i in 1:nrow(colData(sec1_anterior))){
  sample_ID_Sec1Ant[i,1]<-c("Sec1_Anterior")
}
sample_ID_Sec2Post<-matrix(0,nrow(colData(sec2_posterior)),1)
for(i in 1:nrow(colData(sec2_posterior))){
  sample_ID_Sec2Post[i,1]<-c("Sec2_Posterior")
}
sample_ID_Sec1Post<-matrix(0,nrow(colData(sec1_posterior)),1)
for(i in 1:nrow(colData(sec1_posterior))){
  sample_ID_Sec1Post[i,1]<-c("Sec1_Posterior")
}

sample_ID<-rbind(sample_ID_Sec2Ant,sample_ID_Sec1Ant,sample_ID_Sec2Post,sample_ID_Sec1Post)

meta<-cbind(tempcom,sample_ID)

transcrip_edge_num <- 0 
transcrip_adjacency <- list("Node1","Node2")
for(i in 1:(nrow(meta)-1)){
  for(j in (i+1):nrow(meta)){
    if((meta[i,1]==meta[j,1]) & (meta[i,2]!=meta[j,2])){
      transcrip_edge_num <- transcrip_edge_num + 1 
      transcrip_adjacency$Node1[transcrip_edge_num]<-i-1
      transcrip_adjacency$Node2[transcrip_edge_num]<-j-1
    }
  }
}

final_transcrip_adjacency<-cbind(transcrip_adjacency$Node1,transcrip_adjacency$Node2)

```

# Spots Spatial Graph
After extracting the 2D spatial positions of spots for each section and adding the offsets as explained in the MUSTANG manuscript one can use the `get_edges` function from `BayesTME` package to extract the spots spatial graph of the multi-samples ST data using below codes in Python:

```python
from bayestme import data
import pandas
import numpy
from bayestme import utils

temp_pos = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/tot_pos.csv',header=0,index_col=0,sep=",")
pos = temp_pos.to_numpy()

neighbors = utils.get_edges(pos=pos, layout=1)

pandas.DataFrame(neighbors).to_csv("~/Spatial/res_mouseBrain_total/spatial_neighbors.csv",index=True)

```

# Spots Similarity Graph
Now that we have extracted both spots transcriptional and spatial edges we can construct the spots similarity graph edges in R by running below code:

```r
transcrip_neighbors<-read.csv(file = "~/Spatial/res_mouseBrain_total/transcrip_neighbors.csv")
spatial_neighbors<-read.csv(file = "~/Spatial/res_mouseBrain_total/spatial_neighbors.csv",row.names = 1)

colnames(transcrip_neighbors)<-colnames(spatial_neighbors)
similarity_neighbors<-rbind(transcrip_neighbors,spatial_neighbors)

write.csv(similarity_neighbors,file="~/Spatial/res_mouseBrain_total/similarity_neighbors.csv",quote = F,row.names = T)


```
# Bayesian Deconvolution Analysis
Now that we have our all three necessary inputs (i.e. Concatenated spots*genes gene expression matrix, concatenated spots positions with offsets and spots similarity graph), we move to the last step of analyzing our multi-sample ST data with MUSTANG that is performing Bayesian deconvolution analysis. We can do this in Python by running below line of codes:

```python
import pandas
import numpy

temp_counts = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/tot_ST.csv',header=0,index_col=0,sep=",")
counts = temp_counts.to_numpy()

temp_pos = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/tot_pos.csv',header=0,index_col=0,sep=",")
pos = temp_pos.to_numpy()

temp_genes = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/Gene_names.csv',header=0,sep=",")
genes = temp_genes.to_numpy()
genes = genes.reshape(-1)

temp_mask = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/tissue_mask.csv',header=0,sep=",")
mask = temp_mask.to_numpy()
mask = mask.reshape(-1)

temp_barcodes = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/Barcodes.csv',header=0,sep=",")
barcodes = temp_barcodes.to_numpy()
barcodes = barcodes.reshape(-1)

neighbors_tot = pandas.read_csv('~/Spatial/Data/MouseBrain/Total/similarity_neighbors.csv',header=0,index_col=0)
neighbors_tot.reset_index(drop=True, inplace=True)
neigh = neighbors_tot.to_numpy()


stdata= data.SpatialExpressionDataset.from_arrays(raw_counts=counts,
         positions=pos,
         tissue_mask=mask,
         gene_names=genes,
         layout=data.Layout.HEX,
         barcodes=barcodes)

#len(stdata.gene_names)

stdata_stddev_filtered = gene_filtering.select_top_genes_by_standard_deviation(
    stdata, n_gene=1000)

#len(stdata_stddev_filtered.gene_names)

spot_threshold_filtered = gene_filtering.filter_genes_by_spot_threshold(
    stdata_stddev_filtered, spot_threshold=0.95)

stdata_filtered = gene_filtering.filter_ribosome_genes(spot_threshold_filtered)

print('{}/{} genes selected'.format(len(stdata_filtered.gene_names), len(stdata.gene_names)))

from pathlib import Path

best_lambda = 1000
best_n_components = 11

deconvolution_result = deconvolution.deconvolve(
        reads=stdata_filtered.reads,
        edges=stdata_filtered.edges,
        n_gene=1000,
        n_components=best_n_components,
        lam2=best_lambda,
        n_samples=100,
        n_burnin=1000,   # 2000 in publication
        n_thin=5,        # 5 in  publication
        bkg=False,
        lda=False,
        similarity_graph=neigh)


Path("./deconvolution_plots_transcrip").mkdir(exist_ok=True)
Path("./deconvolution_res_transcrip").mkdir(exist_ok=True)

deconvolution.add_deconvolution_results_to_dataset(stdata=stdata_filtered, result=deconvolution_result)

deconvolution_result.save('./deconvolution_res_transcrip/deconvolve_result.h5')

deconvolution.plot_cell_num_scatterpie(
        stdata= stdata_filtered,
        output_path= './deconvolution_plots_transcrip/Scatter_piechart.png')

deconvolution.plot_cell_prob(
        stdata=stdata_filtered,
#        result=deconvolution_result,
        output_dir='./deconvolution_plots_transcrip',
        output_format= 'png',
)

deconvolution.plot_cell_num(
        stdata=stdata_filtered,
#        result=deconvolution_result,
        output_dir='./deconvolution_plots_transcrip',
        output_format= 'png',
)



```
