# Enrichment analysis of SysMalVac results

This folder structure contains all the scripts and data used to produce the enrichment analysis of Moncunill et al. (2019) manuscript.

* Input data, intermediate data, a raw version of the results (together with all the necessary scripts) can be downloaded from Figshare: [sysmalvac_enrichment_scripts.zip]() file.
* Final results can be downloaded from Figshare: [sysmalvac_deliverable.zip]() file.
* An even more succint version of the data is provided in the manuscript as Supplementary Data in the manuscript.

For specific questions about the scripts, please send an email to <miquel.duran@irbbarcelona.org>.
For general questions about the SysMalVac study, please contact <gemma.moncunill@isglobal.org> or <carlota.dobano@isglobal.org>.

## Specifications

This analysis contains many enrichment analyses. Most processes were run with a local SGE cluster at [IRB Barcelona](http://irbbarcelona.org). The main script to set up jobs is `./setupArrayJob.py`. We provide a [Singularity](https://www.sylabs.io/docs/) image containing the necessary dependencies. The folder contains Python (`*.py`), IPython Notebooks (`*.ipynb`) and R-package (`*.R`) scripts.

## Guide through the scripts

1. The `Parsing.ipynb` notebook was used to provide some gene mappings, given the differential gene expression analysis. Input data for the scripts below is provided as `./data`.
2. SGE scripts are denoted by a `_` before the script name. At IRB, queueing systems are called `all.q` and `fast.q`. This should be edited if another cluster is used.
   * `_enrichment.py` handles the [GSEA](http://software.broadinstitute.org/gsea/index.jsp) analysis against modules and gene sets.
   * `_enrichr.py` handles the [EnrichR](https://amp.pharm.mssm.edu/Enrichr/) analysis.
   * `_enrichment_wgcna.py` handles the [GSEA](http://software.broadinstitute.org/gsea/index.jsp) analysis based on co-expression (WGCNA) analysis (see below).
3. Locally-run `R` scripts. These are light-computation scripts that do not require SGE computation.
   * `wgcna.R` performs WGCNA analysis based on the [DCGL](https://cran.r-project.org/web/packages/DCGL/index.html) package.
   * `camera_and_tmod.R` performs Camera and TMod analyses to complement GSEA.
4. Running the scripts above produces a complex folder structure.
   * `sysmalvac_results` folders contain the GSEA, Tmod, Camera and EnrichR analysis of differential expression (conventional enrichment).
   * `sysmalvac_results_wgcna` folders contain the GSEA results of the differential co-expression analysis.
5. Summarizing results.
   * `summarize.ipynb` produces a summary table of the conventional enrichment results.
   * `summarize-wgcna.ipynb` produces a summary table of the WGCNA enrichment results.
6. For convenience, we manually assembled a deliverable (`sysmalvac_deliverable.zip`) containing the relevant input data, conventional enrichment results (`sysmalvac_diff_expr`) and WGCNA enrichment results (`sysmalvac_wgcna`).



