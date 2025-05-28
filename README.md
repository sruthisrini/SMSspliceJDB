# SMSspliceJDB
SMSspliceJDB is a comprehensive and user-friendly database designed to support researchers studying synonymous and non-synonymous mutations with potential splicing implications in cancer patients.

SMSspliceJDB provides a centralized repository of mutations located within acceptor and donor splice junctions, helping researchers understand how these variants may impact RNA splicing and protein function. It includes:

- Experimentally validated mutations
- Predicted splice site variants
- Contextual protein plots
- Downloadable data for in-depth analysis

## Data Inputs
- genes_from_gencode.txt — Gene coordinates from the GENCODE consortium (chromosome, start, end, strand, transcript, gene name).
- variant.txt — Variant data from COSMIC related to six splice site positions (3 acceptor, 3 donor).
- domain_new.txt — Protein domain information associated with genes.
- exp_val_mutations.txt — Experimentally validated mutations associated with splice site positions.

## Features
- Uses GenomicRanges to represent mutations and protein domains as genomic features.
- Visualizes mutations with lollipop plots via the trackViewer package.
- Enables exporting filtered mutation data as CSV files.
- Employs a clean UI theme using shinythemes for a better user experience.

## How to Use
- Enter a gene symbol (e.g., TP53, BRCA1) in the input box on either tab.
- Click "Submit" to generate plots and tables based on available data for the gene.
- Explore the lollipop plots for mutation positions, protein domain context, and mutation distribution.
- View or download the filtered mutation data table for further analysis.

## Visualizations
Lollipop plots display mutation locations and protein structures within user-defined genes, focusing on splice site junctions in cancer patients and experimentally validated splice site mutations.

![image](https://github.com/user-attachments/assets/08c202a0-91d6-4686-99c5-e669be7f5710)

The lollipop plot shown here, AI-generated using DALL·E in ChatGPT, visually represents experimentally validated splice site mutations. The original plot cannot be uploaded due to data privacy and publication restrictions.

## Dependencies
- shiny
- GenomicRanges
- ggplot2
- trackViewer
- maftools
- shinythemes


