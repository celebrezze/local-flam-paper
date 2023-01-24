# local-flam-paper

## A GitHub Repository for: 

## Live fuel moisture and shoot water potential exhibit contrasting relationships with leaf-level flammability thresholds during laboratory flammability tests

## Indra Boving, Joe Celebrezze, Aaron Ramirez, Ryan Salladay, Leander Love-Anderegg and Max Moritz

--------------------------------

## Introductory Statement
This repository is meant for the storage and sharing of data, scripts, figures, and mixed effects model result tables related to the paper referred to above. This paper looks at the relationship between tissue flammability of two dominant chapparal species -- *Ceanothus megacarpus* (Bigpod ceanothus) and *Adenostoma fasciculatum* (chamise) -- and water status (water potential) and content (LFM). In looking at this relationship, it also looks at the effects of different species and sampling years/months on this relationship. Tissue-level flammability was measured using an epiradiator-style chamber over the course of 4 years (2016-2020).

--------------------------------

## Breakdown of Folders

**Raw Data**:
The *raw-data* folder consists of three datasets in .csv format. Below are brief descriptions of each of these dataframes and their use(s); however, for a more comprehensive breakdown of the data, see the metadata.

  *flam.local.alldates.csv*: this includes the flammability testing results for chapparal shrubs, *Adenostoma fasciculatum* and *Ceanothus megacarpus*, reporting a variety of metrics including the sample, round and siteneeded to identify the invidual sample and flammability test (see metadata for exactly what this means), live fuel moisture (LFM), the components necessary to calculate LFM (dry weight and fresh weight), and water potential data along the benchtop drydown; a variety of flammability metrics including flame height (fh), time to first glow (ttfg), glow to ignition (gti), time to ignition (tti), flame duration (fd), glow duration (gd), post-flame glow (pfg), maximum temperature (temp.max), and temperature at ignition (ignition.temp); sample weight; the proportion of new growth, and a variety of other variables. 
  
  *local_field_data.csv*: this includes field-collected values of live fuel moisture (measured at midday) and water potential (measured pre-dawn and at midday). This was only done on one occassion for this study, on 9/18/20, thus it is not a large component of the study and was not discussed in the main text or supplementary index.
  
  *Precip_SBBG.csv*: this includes daily precipitation values as measured at the Santa Barbara Botanical Gardens for the County of Santa Barbara. This data was accessed at: https://www.countyofsb.org/2256/Historical-Rainfall-Reservoir-Informatio .

**Processed Data**:
Versions of the above datasets after wrangling data on *1.0_data_wrangling.Rmd* and *1.1_SBBG_precip.Rmd*. For a comprehensive breakdown of each dataset, see the metadata

**Scripts**:
This folder includes scripts that we used to wrangle data, complete analyses, and design tables and figures for the main text, the supplementary index, and for exploratory analyses. The scripts are grouped by numbers in a logical order which follows the order presented in the manuscript.

  *1.0_data_wrangling.Rmd*: this script involves necessary data wrangling ffor most of our analyses in the following scripts. It works off of the *flam.local.alldates.csv* raw data above and cleans up the data, removing NA values and outliers when necessary, filtering out hot plate burns (this method was used in other studies, but this study focused on epiradiator burns); creates new columns used for different analyses; and designs datasets specific to certain analyses -- the segmented regression analyses (in *4.x*_...Rmd scripts) and the mixed effects models analyses (in *3.x*_...Rmd scripts)
  
  *1.1_SBBG_precip.Rmd*: this script is a data wrangling script specific to the precipitation data in *Precip_SBBG.csv* and extracts the precipitation for two and four months prior to testing by summing up daily precipitation values and then making the dataframe in *processed-data*, *precip.clean.csv*; note: the information gathered from this script was used in *1.0_data_wrangling.Rmd* to add a precip_2mo column (a column describing precipiation (inches) for 2 months prior to sampling month, see metadata for more information); however, the dataframe made in this script was not utilized -- instead it was just the precip_2mo values manually inputted
  
  *2_PCA.Rmd*: this script contains code necessary to run a principal component analysis (PCA) on the flammability data, visualize the PCA results, and create a table of the PCA results. Also, it contains the necessary code for supplementary PCA visualizations, looking at the PCA for each species side-by-side
  
  *3.0_mixed_effects_models.Rmd*: this script, unsurprisingly, showcases code necessary for the linear mixed effects models and associated visualizations and tables. We tried a variety of data manipulations and combinations of variables in the mixed effects models, so the script is rather long, but primary model tables are made obvious in the script with all capital letters "MAIN TABLE" or "MAIN MODEL SELECTION"
  
  *3.1_season_and_species_effects.Rmd*: this script includes code necessary to investigate interspecific and seasonal differences in tissue-level flammability. Like *3.0_mixed_effects_models.Rmd*, this script is rather long and convoluted in some ways; however we tried to make primary visualizations, tables and analyses and supplementary tables or visualizations obvious with headers and included a relatively well hashed-out table of contents associating with headers.
  
**Source Code**:


**Figures**:
