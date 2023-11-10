# Bayesian LCLMM monotonic me
Code accompanying article "A latent class linear mixed model for monotonic continuous processes measured with error" by Espin-Garcia (1), Naranjo (2), and Fuentes-Garcia (2), 2023.

(1) Department of Epidemiology and Biostatistics, University of Western Ontario, Dalla Lana School of Public Health and Department of Statistical Sciences, University of Toronto, and Department of Biostatistics and Schroeder Arthritis Institute, University Health Network, Canada
(2) Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México, Ciudad de México, Mexico.

This repository contains five files and one folder:

### JAGS models
The JAGS code for the described models is located in files:
* `lclmm_Kclasses_lclmm.jag` - latent class linear mixed model pre-label switching
* `lclmm_Kclasses_lclmm_labelswitching.jag` - latent class linear mixed model post-label switching
* `lclmm_Kclasses_error_NOincrease.jag` - proposed model pre-label switching
* `lclmm_Kclasses_error_NOincrease_labelswitching.jag` - proposed model post-label switching

### Simulation code
File "lclmm_Kclasses_monotonic_me_sim.R" produces an analogous analysis to the one described in Section 5 ‘Simulation example’, but with different simulated data and a single set of 20 replicates. 

### AppOAI
This folder contains detailed R scrips to reproduce the analysis performed on the motivating data from the Osteoarthritis Initiative (OAI). Note that no data from the OAI were made available. Researchers can request access to these data through [this link](https://nda.nih.gov/oai/), which containts more information of the study.


###########################################################################

To run the files in a local computer, do the following.
 
1.- Download JAGS from www.mcmc-jags.sourceforge.net/

2.- Install the necessary packages to run the R file as indicated within each file. 

3.- Change the address indicated in `setwd("HERE")`, where "HERE" corresponds to the location of the JAGS codes.

4.- Run the R files within each folder to simulate/process the data, estimate the parameters of the model, and visualize the results.

###########################################################################
