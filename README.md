# sRNA Homolog Heatmap Generator

This project provides a workflow for visualizing the phylogenetic distribution of small RNA (sRNA) homologs. It consists of an R script to prepare the data and a Shiny web application to generate interactive heatmaps.

## Project Overview

The goal of this project is to take sRNA homolog data from multiple FASTA files produced by the GLASSgo tool, process this data to include taxonomic information, and then visualize the distribution of these homologs across different taxonomic levels in a heatmap.

### Key Features

* **Data Preparation**: The `prepare_data.R` script automates the process of extracting information from FASTA files and enriching it with taxonomic data.
* **Interactive Visualization**: The Shiny app (`app.R`) provides a user-friendly interface to filter the data and customize the heatmap appearance.
* **Exportable Results**: Both the generated heatmap (as a PDF) and the underlying data matrix (as a CSV) can be downloaded from the app.


## Data Preparation (`prepare_data.R`)

This script is the first step in the workflow. It needs to be run once to process your raw data and create the `heatdata.Rdata` file required by the Shiny app.

### Prerequisites

* **R**: Make sure you have R installed.
* **taxonkit**: This is a command-line tool required for fetching taxonomy data. It needs to be installed and accessible in your system's PATH. You can find installation instructions [here](https://bioinf.shenwei.me/taxonkit/usage/).
* **GLASSgo fasta files**: The fasta headers need to contain all required information e.g. **taxID**, **percent identitity** to reference sRNA or **genome accession** as provided by GLASSgo.
    * **Example Header:**
        ```
        >CP010281.1:c2139966-2139810 Salmonella enterica subsp. enterica serovar Newport str. CVM 22513, complete genome-p.c.VAL:99.36%-taxID:796732
        ```

### How to Run

1.  Place the `prepare_data.R` script in the root of your project directory.
2.  Ensure all your FASTA files (.fa, or .fasta) are available in the project directory.
3.  Run the script from your R console:
    ```R
    source("prepare_data.R")
    ```
4.  The script will process the files and generate `heatdata.Rdata`.

## Shiny App (`app.R`)

The Shiny app provides an interactive interface for exploring the processed data.

### Prerequisites

* **R and Shiny**: You need R and the `shiny` package installed.
* **Required R Packages**: Install the necessary packages by running this command in your R console:
    ```R
    install.packages(c("shiny", "shinydashboard", "ComplexHeatmap", "colourpicker", "circlize"))
    ```

### How to Run

1.  Make sure the `heatdata.Rdata` file is in the same directory as `app.R`.
2.  Open `app.R` in RStudio and click "Run App", or run the following command in your R console:
    ```R
    shiny::runApp("app.R")
    ```
This will launch the web application in your default browser.

## Using the `example.fa` file

The `example.fa` file is a sample of the expected input format for the fasta files. It contains FASTA entries with headers that include the necessary metadata for the `prepare_data.R` script to parse. You can use this file to understand the required format for your own data.

## References

* **GLASSgo:** Lott, S. C., et al. (2018). "GLASSgo: a fast and flexible algorithm for finding common substructures in a set of RNA 2D structures." *Nucleic Acids Research*. [doi: 10.1093/nar/gkx1268](https://doi.org/10.1093/nar/gkx1268)
* **TaxonKit:** Shen, W., & Xiong, J. (2019). "TaxonKit: a practical and efficient NCBI taxonomy toolkit." *Journal of Genetics and Genomics*. [doi: 10.1016/j.jgg.2019.04.003](https://doi.org/10.1016/j.jgg.2019.04.003)
