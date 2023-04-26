# local-flam-paper

## A GitHub Repository for: 

## Live fuel moisture and shoot water potential exhibit contrasting relationships with leaf-level flammability thresholds during laboratory flammability tests

### Indra Boving, Joe Celebrezze, Aaron Ramirez, Ryan Salladay, Leander Love-Anderegg and Max Moritz

--------------------------------

## Introductory Statement
This repository is meant for the storage and sharing of data, scripts, figures, and linear mixed effects model results tables related to the above paper. This paper looks at the relationship between tissue flammability of two dominant chaparral species -- *Ceanothus megacarpus* (bigpod Ceanothus) and *Adenostoma fasciculatum* (chamise) -- and water status (water potential) and content (live fuel moisture). While looking at this relationship, it also investigates the effects of different species and sampling years/months on this relationship. Tissue-level flammability was measured using an epiradiator-style chamber over the course of 4 years (2016-2020).

--------------------------------

## Table of Contents

[Breakdown of Folders](https://github.com/celebrezze/local-flam-paper#breakdown-of-folders)

[- Raw Data](https://github.com/celebrezze/local-flam-paper#raw-data)
  
[- Processed Data](https://github.com/celebrezze/local-flam-paper#processed-data)
  
[- Scripts](https://github.com/celebrezze/local-flam-paper#scripts)
  
[- Figures](https://github.com/celebrezze/local-flam-paper#figures)
  
[Contact Information](https://github.com/celebrezze/local-flam-paper#contact-information)

--------------------------------

## Breakdown of Folders

### Raw Data:
The **raw-data** folder consists of three datasets in .csv format. Below are brief descriptions of each of these datasets and their use(s); however, for a more comprehensive breakdown of the data, see the metadata.

  *flam.local.alldates.csv*: this includes the flammability testing results for chaparral shrubs, *Adenostoma fasciculatum* and *Ceanothus megacarpus*, reporting a variety of metrics including the sample, round and site needed to identify the individual sample and flammability test (see metadata for exactly what this means), live fuel moisture (LFM), the components necessary to calculate LFM (dry weight and fresh weight), and water potential data along the hydration gradient established via a benchtop drydown; a variety of flammability metrics: flame height (fh), time to first glow (ttfg), glow to ignition (gti), time to ignition (tti), flame duration (fd), glow duration (gd), post-flame glow (pfg), maximum temperature (temp.max), and temperature at ignition (ignition.temp); sample weight; the proportion of new growth; and a variety of other variables. 
  
  *local_field_data.csv*: this includes field-collected values of live fuel moisture (measured at midday) and water potential (measured pre-dawn and at midday). This was only done on one occasion for this study, on 9/18/20, thus it is not a large component of the study and was not discussed in the main text nor the supplementary index.
  
  *Precip_SBBG.csv*: this includes daily precipitation values as measured at the Santa Barbara Botanical Gardens for the County of Santa Barbara. This data was accessed at: https://www.countyofsb.org/2256/Historical-Rainfall-Reservoir-Informatio and wrangled to our needs in *1.1_SBBG_precip.Rmd*

### Processed Data:
This folder includes versions of the above datasets after data wrangling in *1.0_data_wrangling.Rmd* and *1.1_SBBG_precip.Rmd*. It also includes a dataset wrangled prior to the creation of this repository which includes the data garnered from running pressure-volume (PV) curves in *pv_summary_df_timing.csv*. For a comprehensive breakdown of each dataset and descriptions of each variable, see the metadata.

### Scripts:
This folder includes scripts that we used to wrangle data, complete analyses, and design tables and figures for the main text, the supplementary index, and for exploratory analyses. The scripts are grouped by numbers in a logical order which follows the order presented in the manuscript.

  *1.0_data_wrangling.Rmd*: this script involves necessary data wrangling for most of our analyses in the following scripts. It works off of the *flam.local.alldates.csv* raw data above and cleans up the data, removing NA values and outliers when necessary, filtering out hot plate burns (this method was used in other studies and at-times burned samples alongside the epiradiator method, but this study focused solely on epiradiator burns); creates new columns used for different analyses; and designs datasets specific to certain analyses -- the mixed effects models analyses (in *3.x*_...Rmd scripts, datasets = *mem.data.local.epi.alldates.csv* and *mem.data.subset.alldates.csv*) and the segmented regression analyses (in *4.x*_...Rmd scripts, dataset = *seg.reg.data.local.csv*).
  
  *1.1_SBBG_precip.Rmd*: this script is a data wrangling script specific to the precipitation data in *Precip_SBBG.csv* and extracts the precipitation for two and four months prior to testing by summing up daily precipitation values and then saving the dataset in *processed-data*, *precip.clean.csv*. Note: the information gathered from this script was used in *1.0_data_wrangling.Rmd* to add a precip_2mo column (a column describing precipitation (inches) for 2 months prior to sampling month, see metadata for more information); however, the dataset made in this script was not utilized. Instead we manually inputted precip_2mo values, so *1.0_data_wrangling.Rmd* can be ran without issue prior to *1.1_SBBG_precip.Rmd*.
  
  *2_PCA.Rmd*: this script contains code necessary to run a principal component analysis (PCA) on the flammability data, visualize the PCA results, and create a table of the PCA results. Furthermore, it contains the necessary code for supplementary PCA visualizations, looking at the PCA for each species side-by-side. Lastly, it involves some data wrangling to input the resulting principal component values into necessary datasets, as we utilized PC1 and PC2 as proxies of ignitability-combustibility and consumability (respectively) extensively in the analyses in order to make the analysis less convoluted and collapse the multiple flammability metrics into relatively easily understood flammability proxies.
  
  *3.0_mixed_effects_models.Rmd*: this script, unsurprisingly, showcases code necessary for the linear mixed effects models and associated visualizations and tables. We tried a variety of data manipulations and combinations of variables in the mixed effects models, so the script is rather long, but primary model tables are made obvious in the script with all capital letters "MAIN TABLE" or "MAIN MODEL SELECTION". Many of the output tables from this script are saved in *MEM.figures*.
  
  *3.1_season_and_species_effects.Rmd*: this script includes code necessary to investigate interspecific and seasonal differences in tissue-level flammability. Like *3.0_mixed_effects_models.Rmd*, this script is rather long and convoluted in some ways; however we tried to make primary visualizations, tables and analyses and supplementary tables or visualizations obvious with headers and included a relatively well hashed-out table of contents associating with headers. These headers are -- in some cases -- preceded by numbers which correspond with the following processes/analyses: 1. data wrangling, 2. season effects, and 3. species effects.
  
  *4.0_segmented_rand_reg.Rmd*: this script includes code for the segmented regression with random effects analyses (the main segmented regression analyses for this study). As noted in the methods section of the paper, package 'nlme' was used in this script rather than the package 'lme4' used in the *3.0...* and *3.1...* scripts due to its improved compatibility with package 'segmented'. This script has all segmented regression analyses for *C. megacarpus* and *A. fasciculatum* using both water potential and LFM as primary predictors and also has two supplementary analyses using (1) -1/water potential, a transformation of water potential commonly used when conducting pressure-volume curves, as the primary predictor and (2) splitting the dataset by sampling date prior to running segmented regressions. We wrote a function (and modified it slightly for different uses) to visualize the results of the segmented regressions. We had to do this since we used segmented.default instead of segmented.lme when running our segmented regressions due to segmented.lme still being in early stages and plot.segmented.lme producing undesirable results (or no results for us, rather, as our 'individual' variable had too many levels as a random effect in the linear mixed effects models for segmented.lme to run properly). IMPORTANT NOTE: Because of the processing requirements of running these analyses, this script can take especially long to run and often stalls out RStudio for a while when you open it up (in the past, it has taken up to twenty minutes for RStudio to be functional after opening the script depending on the processing power of the computer). That being said, after RStudio levels out, the script runs without errors.
  
  *4.1_segmented_non_random_reg.Rmd*: this script was written before we figured out how to include random effects for the segmented regressions. It contains segmented regression analyses using linear models and serves as somewhat of an exploratory analysis. However, due to it informing our conclusions regarding how we interpreted the results of the segmented regressions with random effects, we felt it would be prudent to include it in the main scripts rather than put it into the *extra-analyses* folder. It is organized differently than *4.0*, as it separates the script by predictor (LFM and water potential) rather than by species and predictor.
  
  *4.2_bootstrapped_seg_reg.Rmd*: this script serves as a supplementary analysis to the segmented regression with random effects. We used bootstrapped data to coerce an even data distribution along the x-axis (LFM or water potential) and included the average threshold (after 100 iterations of segmented regressions per predictor (n = 2), flam. metric (n = 11) and species (n = 2)) in supplementary tables. This was done without random effects, as including random effects caused the script to return errors too frequently to do 100 iterations for any of the flam. metrics (because any time an analysis does not return a threshold, the for loop breaks and causes an error to pop up); however, we still believe this was a useful analysis to conduct to test whether or not the thresholds identified in *4.0* were a result of uneven data distribution across the x-axis. Note: For cases that for loops broke down too frequently due to a threshold not being found on certain iterations, we set up a series of for loops in which we selected seeds that did not cause the for loop to break down. Otherwise, we selected a series of seeds (1:100 for CEME and 151:250 for ADFA that led to the code running smoothly and allowed for our results to be repeatable).
  
  *5.0_comparing_MEM_with_seg_reg.Rmd*: this script compares R^2 values, AIC scores, and BIC scores for the mixed effects models and the segmented regressions with random effects to get a sense of which one predicts the data better. This was used to inform our interpretations of the segmented regression results.
  
  *6.0_LDMC_and_RWC.Rmd*: this script contains code to calculate relative water content (RWC) and leaf dry matter content (LDMC) from the branchlets used for flammability benchtop drydowns using extrapolated saturated water contents. LDMC is compared to TLP to see if these are correlated, RWC is modeled as a predictor of flammability (for both linear and segmented mix-effects regressions), as is LDMC. 
  
  *extra-analyses*: this folder contains scripts for exploratory or extra analyses that we deemed unnecessary for the interpretation of the main text or the supplementary index. They are numbered similarly to the above scripts and belong in those groups. Because these are *extra* analyses, the descriptions of each script as follows will be abridged.
  
  *3.2_pre_post_tlp_MEM.Rmd*: this script contains an extra analysis looking at the results of linear mixed effects models if the data was split by the turgor loss point (TLP) prior to running the linear mixed effects models.
  
  *4.3_segmented_rand_reg_limtedLFM.Rmd*: this script contains an extra analysis running segmented regressions with random effects on a dataset that excludes larger or smaller LFM values than seen in landscape-level monitoring at the Santa Barbara Botanical Gardens
  
  *4.4_segmented_rand_reg_percentile.Rmd*: this script contains an extra analysis running segmented regressions with random effects looking at the mean values of the flammability metric per 10% bin of LFM. This was discarded into the *extra-analyses* rather than used to interpret the segmented regression results, as there were issues with unequal data distribution and some bins only had one or two points while others had more than ten, making the conclusions from this script convoluted and nearly impossible to interpret.
  
  *6.0_extra_figures.Rmd*: finally, this script contains a modge podge of extra figures that didn't end up making the main text or supplementary index (with the notable exception of the PV curve summary table (Table S3) which is made in this script). It is organized based on visualization/analysis.

### Figures: 
This folder contains all figures included in the main text (**main-figures** folder) and the supplementary index (**supp-figures** folder), as well as mixed effect model results tables (**MEM.figures** folder) and extra figures and tables that were not included in the main text or supplementary index but were a part of exploratory analyses or were different visualizations for main analyses (**extra-figures** folder). The figure labels describe the figure and, for main and supplementary figures, they are explicitly labelled matching the figures in the paper. 

--------------------------------

## Contact Information

This GitHub repository was worked on by Indra Boving and Joe Celebrezze.

Indra Boving*: bovingi@ucsb.edu

Joe Celebrezze: celebrezze@ucsb.edu

**correspondence*