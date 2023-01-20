Author:
Feng Yu University of California, Merced fyu9@ucmerced.edu

Each row of Table S1 is the data for a single IDR.
Each column is the ensemble property/entropic force strength of a single IDR.
Values will be N/A, if the data is not analyzed for that dataset. 

The schema of Table S1 csv:

protein: protein name
sequence: protein sequence
length: protein sequence length
dataset: For plotting purpose, we have 4 individual dataset in the SI.
	DisProt, GS, PUMA, UGDH. Detail is described in the paper.
(Reference for the following features: http://pappulab.github.io/localCIDER/)
kappa: charge residue distribution
FCR: fraction of the charge residue
NCPR: net charge per residue
kappa_IL: Hydrophobic residue distribution	
kappa_QN: Hydrophilic residue distribution	
kappa_YF: Aromatic residue distribution	
## We use m3,m2,m1,buffer,p1,p2,p3 to represent the solution condition change.
## These solution condition will be in the column names. I will use 'sc' to represent it.
frame_number_'sc': simulated frames number for each condition.(Before enhanced sampling)
ensemble property of each row including: 
	Ree_'sc'_value: average end-to-end distance
	Ree_'sc'_std: standard deviation
	Rg_'sc'_value: average radius of gyration
	RRg_'sc'_std: standard deviation
	H-bonds_'sc'_value: average H-bond per residue
	H-bonds_'sc'_std: standard deviation
	Helicity_'sc'_value: average helicity
	Helicity_'sc'_std: standard deviation
	Beta-Sheet_'sc'_value: average beta-sheet propensity per residue
	Beta-Sheet_'sc'_std: standard deviation
	Asphericity_'sc'_value: average asphericity
	Asphericity_'sc'_std: standard deviation
## Entropic force strength is represented in the OmegaT/OmegaU: Accesable state ratio (correlate with the entropic force strength)
## The colunm name is in the following format 'OmegaT/OmegaU_d_{distance}_t_{angle}_'sc''

distance: d value mentioned in the paper [0,0.2,0.5] nm
theta: angle theta mentioned in the paper [1,30,60] degree