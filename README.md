\### Bayesian Penalisation and variable selection in REMs using WTC Police Calls Dataset Tutorial



This repository contains the materials for the \*\*Bayesian Penalisation and variable selection in REMs using WTC Police Calls dataset tutorial\*\*, including the R code, data, and outputs.

The goal is to provide a fully reproducible workflow that allows anyone to replicate the analyses and figures step by step.





\## 📂 Repository Structure

Bayesian Penalisation and variable selection in REMs using WTC Police Calls Dataset-Tutorial/

├── data/ # Input datasets (UUsummerschool.rdata)

├── figures/ # Generated figures

├── outputs/ # Model results and text outputs

├── R/ # Supporting R scripts

├── WTCPoliceCalls\_tutorial.R # Main tutorial script

└── README.md # Project instructions







\## 🔧 Requirements



\- \[R](https://cran.r-project.org/) (≥ 4.0.0 recommended)  

\- \[RStudio](https://posit.co/download/rstudio-desktop/) (optional but recommended)  



\### 📦 Required R packages



Install the required packages (if not already installed):


install.packages(c("dplyr", "ggplot2", "tidyr", "remstats", "remstimate", "remify", "brms", "bayesplot"))



\## ▶️ How to Run the Tutorial  



1\. \*\*Clone this repository\*\* (or download as ZIP):  

&nbsp;   

&nbsp;   git clone https://github.com/Ali9972Rostami/Bayesian-Penalisation-in-REMs-using-WTCPoliceCalls-dataset

&nbsp;  .git

&nbsp;   cd Bayesian-Penalisation-in-REMs-using-WTCPoliceCalls-dataset







2\. \*\*Open the R project\*\*  

Open `WTCPoliceCalls-Tutorial.Rproj` in \*\*RStudio\*\*.



3\. \*\*Run the tutorial script\*\*  

Open `WTCPoliceCalls\_tutorial.R` in \*\*RStudio\*\*.  



Run the script step by step, or source it entirely:  

&nbsp; source("WTCPoliceCalls\_tutorial.R")







4\. \*\*Outputs\*\*  

\- Figures will be saved to the `figures/` folder.  

\- Model results and summaries will be saved in the `outputs/` folder.  



---



\## 📊 Linking Code to Figures  

For clarity, the R script contains comments linking code chunks to the figures/results in the tutorial paper.  



Example:  

Code for Figure 3

ggplot(data, aes(x, y)) + geom\_point()







This helps users know which part of the code produces which figure.  



---



\## 📝 Citation  

If you use this tutorial or data, please cite:  



Rostami, A. (2025). \*Bayesian Penalisation in REMs using WTCPoliceCalls dataset.\*  

GitHub repository: \[https://github.com/Ali9972Rostami/Bayesian-Penalisation-in-REMs-using-WTCPoliceCalls-dataset]











