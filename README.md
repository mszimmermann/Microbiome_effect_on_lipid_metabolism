# Microbiome\_effect\_on\_lipid\_metabolism

Pharmacokinetic Multi-Compartment Modeling of Fatty Acid Metabolism in Mice

The multi-compartment pharmacokinetic model of fatty acid metabolism in the mouse to assess gut microbiota effects on fatty acid metabolism parameters. The model contains 11 compartments (duodenum, jejunum, ileum contents, duodenum, jejunum, ileum tissue, colon contents, serum, liver, iWAT, eWAT (inguinal and epididymal white adipose tissue)).

**Software requirements:**
The code is developed with MatLab2019b SimBiology toolbox and requires MatLab version 2019b or later.

**License:**
The code is distributed under the terms of the GNU General Public License (please read copyright\_and\_license and LICENSE files for details.

**Installation:**
There is no specific installation required to run the code. The main workflow is provided in the file Scripts/workflow\_model\_combined\_2FA\_data.m.
Expected run time of the workflow is withing 1h on a standard laptop.

Folder contents:

* Data: input time-resolved measurements of labeled palmitic acid (D5-FA16 or FA16\_0Mz275) and tripalmitin (D31-FA16 or FA16\_0Mz301) across mouse tissues at time points 0, 1, 2 and 6 hours. The data was measured in three mouse groups: germ-free mice, OMM11 community colonized mice, and SPF-mice.

  * Data\\lipidome\_data folder contains measurements of total fatty acids and labelled fatty acids across different tissues. 
  * Data\\proteome\_data folder contains protein expression data from the liver of the three mouse groups
* Figures: Figures depicting model fits, group-specific and general model comparisons, and parameter estimates
* Models: physiology-based pharmacokinetic model used in the workflow
* Output: model parameter estimates and quality assessments.
* Scripts: main modellig workflow and supplementary scripts.

  * workflow\_analyze\_proteome\_data.m contains scripts to preprocess proteomics table from the Data folder and prepare it for comparison with public datasets. Script requires UNIPROT mapping file that can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.15092709. This script produces tables in the Output folder that are requires as input to the workflow\_compare\_proteome\_expression\_atlas.m script. 
  * workflow\_map\_affymetrix\_ids.R is a utility script to map Affymetrix gene IDs to ensemble IDs used to compare across datasets. 
  * workflow\_combine\_expression\_atlas\_data.m script combines expression atlas datasets that can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.15092709 into tables and saves them in the Output folder. These tables are required for te comparison script workflow\_compare\_proteome\_expression\_atlas.m. 
  * workflow\_compare\_proteome\_expression\_atlas.m compares lists of differentially expressed genes across expression atlas datasets to the proteomics fold changes generated in workflow\_analyze\_proteome\_data.m script. workflow\_combine\_expression\_atlas\_data.m needs to run to generate the necessary files in the Output folder. 
  * workflow\_analyze\_FA\_data.m contains scripts to analyze 13C labelled fatty acid profiles in three mouse groups and prepare data for modelling. 
  * workflow\_analyze\_totalFA\_data.m contains scripts to compare total fatty acid profiles across tissues and mouse groups. 

**Expected output**:
Expected output is provided in the Figres and Output folders (model parameter estimates, quality assessments, and plots depicting the experimental data and model fits).



This model is part of the work by Maria Zimmermann-Kogadeeva in collaboration with Josef Ecker and teams.

