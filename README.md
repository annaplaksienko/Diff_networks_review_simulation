This repository contains the code for the simulation study for comparison of several methods for estimation of differential Gaussian graphical models.

*Generate data* folder has the code to generate all pairs of $G^{(1)}$ and $G^{(2)}$ graphs and corresponding precision and covariance matrices $\mathbf{\Omega}^{(1)}$, $\mathbf{\Omega}^{(2)}$, $\mathbf{\Sigma}^{(1)}$, $\mathbf{\Sigma}^{(2)}$ that we used in our simulation. Obtained graphs and matrices are also provided. Replications of the datasets $\mathbf{X}^{(1)}$, $\mathbf{X}^{(2)}$ are too heavy to upload but the code for their generation is provided and results do not vary a lot depending on the realization.

*Methods* folder contains the code we used to produce power-FDR curve for five methods.

*Plotting* folder has the code we used for plotting, as well as the dataset with all the performance results.