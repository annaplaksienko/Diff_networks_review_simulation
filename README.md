
In the manuscript, we have 3 types of the differential graph: random, scale-free and star, produced from different pairs of $G^{(1)}$ and $G^{(2)}$ graphs. For random and scale-free differential graphs, corresponding $G^{(1)}$ is the same. Start with file *Generate_G1_Omega1_Sigma1_scalefree_random.R* to generate it. Changing the `m` parameter there you will achieve two different $G^{(1)}$ sizes from the manuscript: $\approx 200$ and $\approx 400$.

Then, to construct $G^{(2)}$ in a way that corresponding $G^{diff}$ is **random**, use *Generate_G_Omega_Sigma_random.R*. By varying the `niter` parameter in the `rewire(G1_graph, keeping_degseq(niter))` you achieve two different $G^{diff}$ sizes from the manuscript: $\approx 50$ and $\approx 100$ (rerun the function until you get the necessary size).

To construct $G^{(2)}$ in a way that corresponding $G^{diff}$ is **scale-free**, use *Generate_G_Omega_Sigma_scalefree.R*. For 50-200 and 100-400 settings simply remove one subgraph of $G^{(1)}$ graph and one block of $\mathbf{\Omega}^{(1)}$ matrix, for 100-200 settings remove two subgraphs/blocks. For 50-400 setting removing the whole block will produce 100 differential edges instead of 50, so we delete half of the edges of each vertex, starting from the the vertex with the highest degree in a subgraph. 

To construct $G^{(1)}$ and $G^{(2)}$ in a way that their differential is a **star** graph, use *Generate_G_Omega_Sigma_star.R*. Re-run `barabasi.game` until you get a graph with two hubs with a degree around 50 for both, so that you can remove one or both hubs to achieve necessary settings.

Use *Generate_X.R* to generate 50 pairs of datasets $\mathbf{X}^{(1)}$, $\mathbf{X}^{(2)}$ for each pair of graphs $G^{(1)}$, $G^{(2)}$.

Graphs and covariance matrices for all settings are available in *G1_G2_Gdiff_Sigma1_Sigma2_all.Rdata* file.

