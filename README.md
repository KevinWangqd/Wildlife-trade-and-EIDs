The project includes a Bayesian joint ZINB model implemented by the INLA project (https://www.r-inla.org/) and an extreme dependence analysis. 

Analysis code.R is the main code for the ZINB model, with data_total.csv (for all EIDs), data_zoo.csv, and data_nonzoo.csv being the main data of EIDs, wildlife trade, and other variables. 
scaled_matrix.csv is the computed adjacent matrix for the spatial gravity model (see Supplementary).

IMEX_grav_rw.Rdata is the conducted result for the Analysis code. R, one may load it directly.

Extremal dependence.R includes the basic instructions of EVA, using the processed data in EVT.csv.
