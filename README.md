# Hilbert Space Theeory of Cycles

This repository implements harmonic projection for identifying persistent cyclic structure in directed network flows. Given an edge-flow representation of directed interactions, the method constructs boundary operators on a simplicial complex and applies harmonic projection to separate the observed flow into dissipative and harmonic components. The harmonic component lies in the kernel of the 1-Hodge Laplacian and represents stable cyclic organization that cannot be explained by gradient-like propagation or local triangular circulation.

The proof-of-concept idea is given in

Chung, M.K. and Anand, D.A.  and El-Yaagoubi, A.B., Jung, J.-H., Qiu, A., Ombao, H. 2026, [Causality as a Minimum Energy Principle](https://arxiv.org/pdf/2604.17151), Proceedings of the Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC), pages 1-5, publihsed as arXiv:2604.17151. 

(C) 2026 Moo K. Chung, University of Wisconsin-Madison


## Harmonic Projection Simulation

In the accompanying simulation, directed time series are generated from networks with known cyclic ground truth. Harmonic projection is then used to recover the persistent cycle and is compared against Granger causality, structural equation modeling, and time-lagged correlation. The method is designed for cyclic causal discovery in recurrent network systems, including resting-state fMRI connectivity, where feedback and higher-order interactions are common.

Run main MATLAB live script [SIMULATION-causality.mlx](
https://github.com/laplcebeltrami/harmonic/blob/main/SIMULATION-causality.mlx). The underlying Hodge Laplacian based codes came from our [Hodge Laplacian library](https://github.com/laplcebeltrami/hodge/blob/main). 

Simulation study published in 

Chung, M.K. and Anand, D.A.  and El-Yaagoubi, A.B., Jung, J.-H., Qiu, A., Ombao, H. 2026, [Causality as a Minimum Energy Principle](https://arxiv.org/pdf/2604.17151), Proceedings of the Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC), pages 1-5, publihsed as arXiv:2604.17151. 


