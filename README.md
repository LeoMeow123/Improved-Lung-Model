# Improved Lung Model on Effects of Pressure- and Volume-controlled Ventilation on Pulmonary Ventilation and Gas Exchange in Obstructive Lung Disease
This is the code for a research paper on Improving the gas exchange lung model conducted by Yipeng Li (y9li@ucsd.edu) and Chenyu Li (chl230@ucsd.edu) from UCSD Bioengineering.

Author: Yipeng Li

VCVandPCV_hold.m is the literature lung model from Tianya Liu et al. Once you open it and run you will have the literature lung model.

For the Improved lung model, enter the file `Improved Model` and Click on VCVandPCV_hold.m. If you are a Windows user, you can directly run the code. if you are a Mac or Linux user, you need to change 
`delete('Vc_VCV.txt');` into `rm('Vc_VCV.txt');`
`delete('Q_VCV.txt');` into `rm('Q_VCV.txt');`
`delete('Vc_PCV.txt');` into `rm('Vc_PCV.txt');`
`delete('Q_PCV.txt');` into `rm('Q_PCV.txt');`
The file has 14 figure outputs, the first 13 are silenced in the provided code. You can manually adjust the figures you want to use.

Statistical Analyses are conducted with both MATLAB (R2024a) and R (4.2.3). Kept in statistical_analysis.m and Stat.qmd. Same as the main file, supplemental figures are silenced.
