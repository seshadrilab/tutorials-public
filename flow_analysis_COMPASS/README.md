# Flow Cytometry R Analysis Tutorial: COMPASS

A tutorial for analyzing flow cytometry data in R using COMPASS (Combinatorial Polyfunctionality Analysis of Single Cells). 

# Introduction

This tutorial contains a basic workflow for analyzing flow cytometry data using COMPASS after gating in FlowJo. COMPASS is a statistical framework that allows for unbiased analysis of antigen-specific T-cell subsets. COMPASS uses a Bayesian hierarchial framework to model all obeserved cell-subsets and select the most likely to be antigen-specific while regularizing the small cell counts that often arise in multi-parameter space. The model gives a posterior probability of specificity for each cell-subset and each sample, which can be used to profile a subject's immune response to external stimuli (e.g., infection or vaccination). 

For this tutorial, we will be analyzing intracellular cytokine staining (ICS) data from a convalescent cohort of COVID-19 subjects who were either hospitalized or not hospitalized. For more information about this cohort, please see the accompanying paper published in JCI Insight on Feb 23, 2021: 
* [Comorbid illnesses are associated with altered adaptive immune responses to SARS-CoV-2](https://pubmed.ncbi.nlm.nih.gov/33621211/)  

For more information about COMPASS, please refer to the original manuscript and documentation:
* https://www.nature.com/articles/nbt.3187
* https://bioconductor.org/packages/COMPASS/

# Installation

This pipeline should be completed in R and RStudio. The following R packages are required for this tutorial:
* [here](https://cran.r-project.org/package=here) 
* [CytoML](https://bioconductor.org/packages/CytoML/)
* [flowCore](https://bioconductor.org/packages/flowCore/)
* [flowWorkspace](https://bioconductor.org/packages/flowWorkspace/)
* [COMPASS](https://bioconductor.org/packages/COMPASS/)

To install the R packages, open an R session and enter the following command lines:
```R
install.packages("here")
install.packages("BiocManager")
BiocManager::install("CytoML")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("COMPASS")
```
# Directory Structure
To directly use the code in this tutorial, you should create an RStudio project and set up your project directory as follows. The `data/` folder should contain all `.xml` and `.fcs` files and the `out/` folder will contain the GatingSet and all COMPASS outputs.
![image](https://user-images.githubusercontent.com/89667908/147301852-f5c1d505-cb04-4841-bdbe-981b0d4bc6f9.png)

You can achieve this directory structure by running the following command lines:
```R
if(!dir.exists(here::here("data"))) {
  cat(sprintf("Creating folder %s\n", here::here("data")))
  dir.create(here::here("data"), recursive = T)
}
if(!dir.exists(here::here("out"))) {
  cat(sprintf("Creating folder %s\n", here::here("out")))
  dir.create(here::here("out"), recursive = T)
}
if(!dir.exists(here::here("out/GatingSet"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/GatingSet")))
  dir.create(here::here("out/GatingSet"), recursive = T)
}
if(!dir.exists(here::here("out/COMPASSResult"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/COMPASSResult")))
  dir.create(here::here("out/COMPASSResult"), recursive = T)
}
```

# Workflow Overview
1. Load data
2. Create a GatingSet
3. Create a COMPASSContainer
4. Run COMPASS
5. Visualize

## Load data
```R
# Load libraries into your current R session
library(here)
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(COMPASS)
```
Download the .xml file and folder containing the associated .fcs files from [Dropbox](https://www.dropbox.com/sh/snr4ycge3jdm50k/AAA4RaV64vDgRXr04wkm9YeMa?dl=0). Then, drag them into the "data" folder of the project directory.

This data has been gated in FlowJo v9 and will be parsed using flowWorkspace.

**FYI:** When using FlowJo v9, the FlowJo workspace must be exported as an .xml file to create a flowjo_workspace object with the function open_flowjo_xml(). However, when using FlowJo v10, the FlowJo workspace can be loaded directly as a .wsp file using the same function open_flow_xml().
```R
# Location of XML file
xml_path <- here::here("data/20200605_COVID_ICS-B3-trunc.xml")
# Location of .fcs files
fcs_subfolder <- here::here("data/20200605_COVID_ICS-B3-FCS-trunc/")
```

Create a flowjo_workspace object with the function open_flowjo_xml().
```{r}
ws <- open_flowjo_xml(xml_path)
```

## Create a GatingSet
### Set-up
A GatingSet holds a set of GatingHierarchy objects, representing a set of samples and the gating scheme associated with each.
Look at the workspace metadata to choose which keywords to extract into the GatingSet. The flowjo_to_gatingset() function parses a flowJo Workspace to generate a GatingSet object.
```R
# Look at all of the keywords
names(fj_ws_get_keywords(ws, 117)) 
# Choose which keywords to keep
keywords2import <- c("EXPERIMENT NAME",
                       "$DATE",
                       "SAMPLE ID",
                       "PATIENT ID",
                       "STIM",
                       "WELL ID",
                       "PLATE NAME") 
sampleGroup <- "Samples"
gs <- flowjo_to_gatingset(ws,                                    
                          name = sampleGroup, 
                          keywords = keywords2import,
                          path = fcs_subfolder, 
                          extend_val = -10000)
```

### Quality Control (QC)
Make sure that the gating trees are consistent for all samples.
```R
pop_lists <- lapply(gs, gh_get_pop_paths)
unique(pop_lists)
```
```R
## [[1]]
## [1] "root"                                     
## [2] "/Time"                                    
## [3] "/Time/LD-3+"                              
## [4] "/Time/LD-3+/1419-3+"                      
## [5] "/Time/LD-3+/1419-3+/S"                    
## [6] "/Time/LD-3+/1419-3+/S/Lymph"              
## [7] "/Time/LD-3+/1419-3+/S/Lymph/4+"           
## [8] "/Time/LD-3+/1419-3+/S/Lymph/4+/107a"      
## [9] "/Time/LD-3+/1419-3+/S/Lymph/4+/154"       
## [10] "/Time/LD-3+/1419-3+/S/Lymph/4+/CCR7+"     
## [11] "/Time/LD-3+/1419-3+/S/Lymph/4+/CD45RA+"   
## [12] "/Time/LD-3+/1419-3+/S/Lymph/4+/IFNG"      
## [13] "/Time/LD-3+/1419-3+/S/Lymph/4+/IL2"       
## [14] "/Time/LD-3+/1419-3+/S/Lymph/4+/IL17"      
## [15] "/Time/LD-3+/1419-3+/S/Lymph/4+/IL4513"    
## [16] "/Time/LD-3+/1419-3+/S/Lymph/4+/TNF"       
## [17] "/Time/LD-3+/1419-3+/S/Lymph/8+"           
## [18] "/Time/LD-3+/1419-3+/S/Lymph/8+/IFNG"      
## [19] "/Time/LD-3+/1419-3+/S/Lymph/CD38+"        
## [20] "/Time/LD-3+/1419-3+/S/Lymph/HLADR+"       
## [21] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+"        
## [22] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a"   
## [23] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154"    
## [24] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+"  
## [25] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+"
## [26] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IFNG"   
## [27] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2"    
## [28] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17"   
## [29] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513" 
## [30] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF" 
```
Remove channels from flow data that are not used by gates.
```R
gs <- gs_remove_redundant_channels(gs) 
```
```R
## drop SSC-H
```
Add names to all channels or change their names.
```R
dput(unname(pData(parameters(gh_pop_get_data(gs[[1]])))[,2]))
```
```R
## structure(c(NA, NA, NA, NA, "CD8b BB700", "TNFa FITC", "CD107a PE-Cy7", 
## "CD154 PE-Cy5", "CD3 ECD", "IL2 PE", "CD4 APC-H7", "IL17a Ax700", 
## "IL4/5/13 APC", "CD14/CD19 BV785", "CCR7 BV711", "CD38 BV605", 
## "L/D", "IFNg V450", "CD45RA BUV737", "HLADR BUV395"), class = "AsIs")
```
```R
markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a",
                 "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19",
                 "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames) <- pData(parameters(gh_pop_get_data(gs[[1]])))[,1]
markernames(gs) <- markernames
pData(parameters(gh_pop_get_data(gs[[1]])))[,c(1,2)]
```
```R
##          name     desc
## $P1      Time     Time
## $P2     FSC-A    FSC-A
## $P3     FSC-H    FSC-H
## $P4     SSC-A    SSC-A
## $P6  <B710-A>     CD8b
## $P7  <B515-A>     TNFa
## $P8  <G780-A>   CD107a
## $P9  <G660-A>    CD154
## $P10 <G610-A>      CD3
## $P11 <G575-A>      IL2
## $P12 <R780-A>      CD4
## $P13 <R710-A>    IL17a
## $P14 <R660-A> IL4_5_13
## $P15 <V780-A>  CD14_19
## $P16 <V710-A>     CCR7
## $P17 <V610-A>     CD38
## $P18 <V510-A>       LD
## $P19 <V450-A>     IFNg
## $P20 <U730-A>   CD45RA
## $P21 <U395-A>    HLADR
```
In this tutorial, we will only run COMPASS on CD4+ T cells, but what if we wanted to also run COMPASS on CD8+ T cells? Currently, most of the gates are missing from the "8+" node. Let's grab the missing gates from the "NOT4+" node and add them under the "8+" node.
```R
gates_to_copy <- c("/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154", 
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+")
for(path in gates_to_copy) {
  gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y=path),
             parent = "/Time/LD-3+/1419-3+/S/Lymph/8+")
}
recompute(gs)
```
```R
# done!
```
Next, let's get rid of the "NOT4+" node and all of its descendants.
```R
gs_pop_remove(gs, "/Time/LD-3+/1419-3+/S/Lymph/NOT4+")
```
Plot the gating tree.
```R
plot(gs, fontsize=15, bool=T)
```
![image](https://user-images.githubusercontent.com/89667908/151273845-6a500bf4-9546-4a3d-880a-0fd08d916a6c.png)


### Save GatingSet
**Note:** this can take a while.
```R
save_gs(gs, here::here("out/GatingSet"), overwrite = TRUE)
```

## Create a COMPASSContainer
A COMPASSContainer is the data structure used to hold data from an ICS experiment. The input for this code is a GatingSet or a GatingSetList. Counts, metadata, and 
single cell data are extracted and fed into the COMPASSContainer constructor.
```R
# Set the seed
set.seed(123)
# A regular expression to match a single node in the gating tree
parent_node <- "4+"
# A character identifying the subject id column in the GatingSet metadata
id <- "SAMPLE ID"
markernames(gs) 
```
```R
##   <B710-A>   <B515-A>   <G780-A>   <G660-A>   <G610-A>   <G575-A>   <R780-A>   <R710-A> 
##     "CD8b"     "TNFa"   "CD107a"    "CD154"      "CD3"      "IL2"      "CD4"    "IL17a" 
##   <R660-A>   <V780-A>   <V710-A>   <V610-A>   <V510-A>   <V450-A>   <U730-A>   <U395-A> 
## "IL4_5_13"  "CD14_19"     "CCR7"     "CD38"       "LD"     "IFNg"   "CD45RA"    "HLADR" 
```
```R
# markermap contains the output of markernames(gs)
markermap <- list("IL2", "IL4_5_13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")
# Assign names to the list of markers based on the gate names
names(markermap) <- paste0(parent_node, "/", c("IL2", "IL4513", "IFNG",
                                               "TNF", "IL17", "154", 
                                               "107a"))
```
Construct the COMPASSContainer. If the number of parent cells is less than countFilterThreshold, we drop that file (default is 5000 cells).
```R
CC <- COMPASSContainerFromGatingSet(gs,
                                    node = parent_node,
                                    individual_id = id,
                                    mp = markermap,
                                    countFilterThreshold = 5000)
```
```R
## Extracting cell counts
## Fetching 4+
## Fetching child nodes
## common markers are: 
## Time FSC-A FSC-H SSC-A CD8b TNFa CD107a CD154 CD3 IL2 CD4 IL17a IL4_5_13 CD14_19 CCR7 CD38 LD IFNg CD45RA HLADR 
## Extracting single cell data for 4+/IL2|4+/IL4513|4+/IFNG|4+/TNF|4+/IL17|4+/154|4+/107a
## ..........................................................................................................................................Creating COMPASS Container
## Filtering low counts
## Filtering 0 samples due to low counts
```
Look at some basic info about our COMPASSContainer.
```R
CC
```
```R
## A COMPASSContainer with 30 samples from 5 individuals, containing data across 7 markers.
```

## Run COMPASS
Fit the COMPASS model using the COMPASSContainer. To fit the COMPASS model, we need to specify how to identify the samples that are our treatment condition and our control condition based on the metadata. Here, we will run COMPASS on the samples stimmed by spike 1 with DMSO as our negative control. For now, let's just do 100 iterations for speed.
```R
fit <- COMPASS(CC,
               treatment = STIM == "Spike 1",
               control = STIM == "DMSO",
               iterations = 100)
```
```R
## There are a total of 5 samples from 5 individuals in the 'treatment' group.
## There are a total of 5 samples from 5 individuals in the 'control' group.
## The model will be run on 5 paired samples.
## The category filter has removed 38 of 54 categories.
## There are a total of 16 categories to be tested.
## Initializing parameters...
## Computing initial parameter estimates...
## Keeping 100 iterations. We'll thin every 8 iterations.
## Burnin for 100 iterations...
## Sampling 800 iterations...
## Done!
## Computing the posterior difference in proportions, posterior log ratio...
## Done!
```
Save the COMPASS run output.
```R
saveRDS(fit, file.path(here::here("out/COMPASSResult"), "COMPASSResult.rds"))
```
Save the Functionality and Polyfunctionality Scores.
```R
FS <- FunctionalityScore(fit)
FS_df <- data.frame(tmp = names(FS), FS = FS)
colnames(FS_df) <- c("SAMPLE ID", "FS")
PFS <- PolyfunctionalityScore(fit)
PFS_df <- data.frame(tmp = names(PFS), PFS = PFS)
colnames(PFS_df) <- c("SAMPLE ID", "PFS")
FS_PFS_df <- merge(FS_df, PFS_df, by = "SAMPLE ID")
write.table(FS_PFS_df,
            file = file.path(here::here("out/COMPASSResult"), "FS_PFS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
FS_PFS_df
```
```R
##   SAMPLE ID         FS        PFS
## 1      133C 0.02811024 0.01214286
## 2      142C 0.03078740 0.01504422
## 3      150C 0.03220472 0.01426871
## 4        23 0.04157480 0.01965986
## 5        25 0.03188976 0.01413265
```

## Visualize
Plot a heatmap of the mean probability of response.
```R
plot(fit, show_rownames = TRUE)
```
```R
## The 'threshold' filter has removed 8 categories:
## IL2&!IL4_5_13&!IFNg&!TNFa&!IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&IFNg&!TNFa&!IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&!IFNg&TNFa&!IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&!IFNg&!TNFa&IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&!IFNg&!TNFa&!IL17a&CD154&!CD107a, IL2&IL4_5_13&!IFNg&!TNFa&!IL17a&!CD154&!CD107a, IL2&!IL4_5_13&!IFNg&TNFa&!IL17a&!CD154&!CD107a, !IL2&IL4_5_13&!IFNg&!TNFa&IL17a&!CD154&!CD107a
```
![image](https://user-images.githubusercontent.com/89667908/151413096-d42a8bfe-6319-43ba-9146-2f5a6a1de3f4.png)
