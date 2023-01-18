Author:
Feng Yu University of California, Merced fyu9@ucmerced.edu

The schema of ensemble data csv:

Protein: Protein name
datatype : datatype of each row including: 
	re: average end-to-end distance
	rg: average radius of gyration
	hb: average H-bond per residue
	helix: average helicity
	beta: average beta-sheet propensity per residue
	asphericity: average asphericity
prot_solv_inter: relative protein:solvent interaction strength
value: value of certain datatype
Std: standard deviation of the previous value
frame_counts: Number of frames
Sequence: protein sequence(all following features will have valuewhen datatype=='feature')	
length: protein sequence length
(Reference for the following features: http://pappulab.github.io/localCIDER/)
kappa: charge residue distribution
FCR: fraction of the charge residue
NCPR: net charge per residue
Hydrophobicity:	hydrophobic residue per residue
Expanding: expanding residue per residue
delta: charge residue distribution(raw data)
kappa_IL: Hydrophobic residue distribution	
kappa_QN: Hydrophilic residue distribution	
kappa_YF: Aromatic residue distribution	


The schema of entropic force csv:
Protein: Protein name
distance: d value mentioned in the paper
theta: angle theta mentioned in the paper
prot_solv_inter: relative protein:solvent interaction strength
OmegaT/OmegaU: Accesable state ratio (correlate with the entropic force strength)
frame_counts: Number of frames