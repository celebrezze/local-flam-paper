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
  
  *Precip_SBBG*: this includes daily precipitation values as measured at the Santa Barbara Botanical Gardens for the County of Santa Barbara. This data was accessed at: https://www.countyofsb.org/2256/Historical-Rainfall-Reservoir-Informatio .

**Processed Data**:
Versions of the above datasets after wrangling data on *1.0_data_wrangling.Rmd* and *1.1_SBBG_precip.Rmd*. For a comprehensive breakdown of each dataset, see the metadata

**Scripts**:
This folder includes scripts that we used to wrangle data, complete analyses, and design tables and figures for the main text, the supplementary index, and for exploratory analyses. The scripts are grouped by numbers in a logical order which follows the order presented in the manuscript.

**Source Code**:


**Figures**:
