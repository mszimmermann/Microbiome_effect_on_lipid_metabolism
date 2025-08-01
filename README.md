# Microbiome_effect_on_lipid_metabolism

Pharmacokinetic Multi-Compartment Modeling of Fatty Acid Metabolism in Mice

The multi-compartment pharmacokinetic model of fatty acid metabolism in the mouse to assess gut microbiota effects on fatty acid metabolism parameters. The model contains 11 compartments (duodenum, jejunum, ileum contents, duodenum, jejunum, ileum tissue, colon contents, serum, liver, iWAT, eWAT (inguinal and epididymal white adipose tissue)).

**Software requirements:**
The code is developed with MatLab2019b SimBiology toolbox and requires MatLab version 2019b or later.

**License:** 
The code is distributed under the terms of the GNU General Public License (please read copyright_and_license and LICENSE files for details.

**Installation:**
There is no specific installation required to run the code. The main workflow is provided in the file Scripts/workflow_model_combined_2FA_data.m.
Expected run time of the workflow is withing 1h on a standard laptop. 

Folder contents:
- Data: input time-resolved measurements of labeled palmitic acid (D5-FA16 or FA16_0Mz275) and tripalmitin (D31-FA16 or FA16_0Mz301) across mouse tissues at time points 0, 1, 2 and 6 hours. The data was measured in three mouse groups: germ-free mice, OMM11 community colonized mice, and SPF-mice.
- Figures: Figures depicting model fits, group-specific and general model comparisons, and parameter estimates
- Models: physiology-based pharmacokinetic model used in the workflow 
- Output: model parameter estimates and quality assessments. 
- Scripts: main modellig workflow and supplementary scripts.

**Expected output**:
  Expected output is provided in the Figres and Output folders (model parameter estimates, quality assessments, and plots depicting the experimental data and model fits). 


This model is part of the work by Maria Zimmermann-Kogadeeva in collaboration with Josef Ecker and teams. 

