# Quantification and Analysis of *in situ* hybridization on Human Tissue Microarrays (TMAs)

This project is for the analysis of mRNA quantifciation of a gene in colorectal cancer TMAs.

## Data set up

The TMAs should be stained according to the optimized protocol from Affymetrix's ViewRNA kit.

Currently, the TMAs are set up with each patient having up to three replicate spots. One picture is taken of each spot. Two to three regions of interest are drawn on each picture using the MetaMorph program. Staining is quantified using the Integrated Morphometry Analysis feature within MetaMorph. Details are in the quantification file.

The data should be cleaned up in Excel and formatted into the following columns:

Image | Probe | TMA | Stage | Sample | Spot | Region | Area | Count | Threshold
imageName.tif | No Probe | TMA1 | Normal | NC1 | A | 1 | 58136 | 6 | 195


Save the files in .csv format.

To do an analysis, you will need both the No Probe and Probe of Interest files.


## Analysis

Analysis.R will provide the intial analysis on the data. It will compute the average mRNA per area for each stage and subtract background (No Probe). A table will be compiled and will be saved in the working directory.

It will then display graphs with the distribution of all stages and stages independently.

The script will then combine the data, perform tests for normality, and conduct a One-Way ANOVA.

A graph will be produced comparing all stages.

## Running the analysis

To run the analysis, you must set a working directory in R that directs to your files.
The files should be named noprobe.csv, probe.csv (ex. pgc1beta.csv)
For the time being, you will have to set the probe file manually in the script. Hopefully I can add prompts in the future.
Once that is set, you can just run the entire script.
