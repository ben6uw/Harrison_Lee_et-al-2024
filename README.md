# Harrison_Lee_et-al-2024
Wide ranging genetic variation in sensitivity to rapamycin in Drosophila melanogaster

This GitHub is organized into parent folders for each of the following code sections:
Each contains at least one .R file with the R code, with reference to the data that is also available within the folder.
The figures [Fig_X] resulting from each analysis are indicated.

1. Phenotype_files
file = phenotype_summary.R
a. Compile and characterize data [Figs 1A, S1]
b. Dose-response [Fig 1b]						
		c. Size analysis [Fig 1c, Fig S2]
	file = food deprivation files/food_deprivation_analysis.R
		d. Food deprivation analysis [Fig 4C] 

3. GWAS_files
	file = GWAS_phenotype.R
		a. GWAS phenotype
	file = GBLUP.R
		b. Heritability [Fig 2A]
		c. Covariance Association Test [Fig 2B, 2C]

4. Metabolomics_files
	file = larva_rapa_metabolome.R
		a. Data normalization
		b. Association with phentoype [Figs 3A, 3B]
		c. Enrichment analysis
	file= starving_larva_Jouandin_Science_2022
		a. Starvation analysis [Fig 4A, 4B]
