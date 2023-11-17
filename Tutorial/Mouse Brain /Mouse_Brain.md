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

