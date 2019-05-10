# Enrichment analysis of SysMalVac gene expression measurements

This repository contains all the scripts used to produce the enrichment analysis of Moncunill et al. (2019) manuscript.

* Input data, intermediate data and a raw version of the results can be downloaded from [Figshare](https://figshare.com/account/home#/projects/63467).
* Final enrichment results can be downloaded from Figshare, too: [sysmalvac_deliverable.zip](https://figshare.com/s/cd96f43063f2ab2c3e5d) file.
* An even more succint version of the data is provided in the manuscript as Supplementary Data in the manuscript.

For specific questions about the scripts, please send an email to <miquel.duran@irbbarcelona.org>.
For general questions about the SysMalVac study, please contact <gemma.moncunill@isglobal.org> or <carlota.dobano@isglobal.org>.

## Specifications

All the scripts necessary to run the analysis are provided within this repository. The procedures involve many enrichment analyses. Most processes were run with a local SGE cluster at [IRB Barcelona](http://irbbarcelona.org). The main script to set up jobs is `./setupArrayJob.py`. We provide a [Singularity](https://www.sylabs.io/docs/) image containing the necessary dependencies. The folder contains Python 2.7 (`*.py`), IPython Notebooks (`*.ipynb`) and R-package (`*.R`) scripts. Scripts are prepared to run with Linux systems.

This repository is provided to ensure data transparency. The outcome of the scripts is conveniently provided as compressed files. To complete the full repository in your home directories, please download the corresponding files from Figshare and uncompress them inside the cloned repository folder:
* [data.tar.gz]()
* [sysmalvac_results.tar.gz]()
* [sysmalvac_results_wgcna.tar.gz]()

## Guide through the scripts

1. The `Parsing.ipynb` notebook was used to provide some gene mappings, given the differential gene expression analysis. Input data for the scripts below is provided as `./data`.
2. SGE scripts are denoted by a `_` before the script name. At IRB, queueing systems are called `all.q` and `fast.q`. This should be edited at the end of the scripts if another cluster is used.
   * `_enrichment.py` handles the [GSEA](http://software.broadinstitute.org/gsea/index.jsp) analysis against modules and gene sets.
   * `_enrichr.py` handles the [EnrichR](https://amp.pharm.mssm.edu/Enrichr/) analysis.
   * `_enrichment_wgcna.py` handles the [GSEA](http://software.broadinstitute.org/gsea/index.jsp) analysis based on co-expression ([WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)) analysis (see below).
3. Locally-run `R` scripts. These are light-computation scripts that do not require SGE computation.
   * `wgcna.R` performs WGCNA analysis based on the [DCGL](https://cran.r-project.org/web/packages/DCGL/index.html) package.
   * `camera_and_tmod.R` performs Camera and TMod analyses to complement GSEA.
4. Running the scripts above produces a complex folder structure.
   * `sysmalvac_results` folders contain the GSEA, [Tmod](https://cran.r-project.org/web/packages/tmod/index.html), [Camera](http://bioconductor.org/packages/release/bioc/html/CAMERA.html) and EnrichR analysis of differential expression (conventional enrichment).
   * `sysmalvac_results_wgcna` folders contain the GSEA results of the differential co-expression analysis.
5. Summarizing results.
   * `summarize.ipynb` produces a summary table of the conventional enrichment results.
   * `summarize-wgcna.ipynb` produces a summary table of the WGCNA enrichment results.

For convenience, we manually assembled a deliverable (`sysmalvac_deliverable.zip`) containing the relevant input data, conventional enrichment results (`sysmalvac_diff_expr`) and WGCNA enrichment results (`sysmalvac_wgcna`).
