
# MC-2020
---
title: "README.md"
author: "Sarah I. Murphy, Tim Lott"
date: "01/07/21"
---

## Purpose: 
This repository is for Milk-MC team to work on modeling projects. 

## Background: 
The first published version of FSL-MQIP fluid milk MC model was published here: [Buehler, A. J., Martin, N. H., Boor, K. J., & Wiedmann, M. (2018). Psychrotolerant spore-former growth characterization for the development of a dairy spoilage predictive model. Journal of Dairy Science.](https://doi.org/10.3168/jds.2018-14501). Michael Phillips edited the original fluid milk MC model, incorporating dynamic temperatures that change according to user-defined stages; M. Phillips also created a respository (FSL-MQIP/MC-Spore) in FSL-MQIP's GitHub organization where he uploaded the model files. On 1/28/20, Sarah Murphy created a new repository (FSL-MQIP/MC-2020) where she uploaded a copy of files stored in FSL-MQIP/MC-Spore. 

## Rules for people working on projects in MC-2020
1. Create a folder for your specific project within MC-2020, including the name of the person leading the project and the name of the project (e.g., Tim_ESL). Once this project folder is created, update this main README.md file in the MC-2020 folder in the Contents section.

2. Within the folder you created, you must have (i) the model RMD file (e.g., MCsim_TimESLModel_010721.Rmd), (ii) the "UtilityFunctions.R" file containing the functions used for your model including descriptions of each of these functions, (iii) a folder named "InputFiles", which contains the files that are inputs for your model, (iv) a folder named "Results" where model results will go, including data and plots, (v) an "Archive" folder (this is where older versions of files will go, or unneeded files, DO NOT DELETE FILES just move them here), (vi) and any other folders as necessary. 

3. Importantly, at the top of your model RMD file, you MUST include information on (i) the original source of your model RMD file, including who authored that file and (ii) an overall summary of your model RMD file (this can be revised over time, but must be complete before publication). You also must include at the top of your model RMD file, "Title", "Author", and "Date", where "Date" must be updated every time you revise the file. 

4. Within the model RMD file, you must include a log of your revisions include that date of your revision, followed by the names of the person or people who made that revision, and a summary of revisions made. You can also include notes for what you want to do for next steps, or issues you are having, etc. as a part of your revisions statement. For each day that you make revisions, when you commit the revisions through GitHub, it's recommended that you copy-paste your revisions statement as the commit comment.

## Contents:
* Files in the main folder of MC-2020 are written so they can be adapted for current and future MC projects. 

* Subfolders in MC-2020 are created for each project. This is the current list of projects within the MC-2020 folder:
- Al_FFAR_SCCSporeModel
- Ariel_2017MC_SporeModel
- Ariel_2018Chobani_YogurtModel
- Forough_FFAR_OptimizationModel
- Sam_FFAR_PPCmodel
- Sarah_CombinedModel
- Sarah_UpstateBuffalo_BactofugeSimulations
- Tim_ESL
- Tim_SchoolMilk
        
    
            




