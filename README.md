# Harmonic Mean Iteratively Reweighted Least Squares for Low-Rank Matrix Recovery

This repository contains MATLAB code to implement a basic variant of the Harmonic Mean Iteratively Reweighted Least Squares (HM-IRLS) algorithm for low-rank matrix recovery, in particular for the low-rank matrix completion problem, and to reproduce the experiments described in the paper:

> C. Kümmerle, J. Sigl.
> "Harmonic Mean Iteratively Reweighted Least Squares for Low-Rank Matrix Recovery", to appear in the Journal of Machine Learning Research (JMLR).
> Available online: https://arxiv.org/abs/1703.05038

The main file is `HM_IRLS.m`. See also the example scripts:
* `script_mc_comparisons.m` - Comparison script between HM-IRLS and two other algorithms on random matrix completion data
* `script_small_example_IRLSvariants.m` - Script illustrating the small example of Section 3 of the paper, comparing HM_IRLS with other IRLS variants for the problem
* `script_HM_IRLS_Figure3.m`  - Script reproducing experiment of Figure 3 of the paper (convergence rates of HM-IRLS and other IRLS variants for easy problems)
* `script_HM_IRLS_Figure4.m`  - Script reproducing experiment of Figure 4 of the paper (convergence rates of HM-IRLS and other IRLS variants for hard problems)
* `script_HM_IRLS_Figure5.m`  - Script reproducing experiment of Figure 5 of the paper (convergence rates of HM-IRLS and other IRLS variants for very hard problems)

## Version history
* Version 1.1, updated 10/01/2018
* Version 1.0, 3/14/2017

## Author
Christian Kümmerle ([website](http://www-m15.ma.tum.de/Allgemeines/ChristianKuemmerle)) 

## Acknowledgments
For the purpose of comparison with other popular algorithmic approaches, we included code by
* Bart Vandereycken ("LRGeomCG"), corresponding to the paper "Low-Rank Matrix Completion by Riemannian Optimization", SIAM J. Optim., 23(2), 1214–1236.
* Zaiwen Wen, Wotao Yin and Yin Zhang ("LMaFit"), corresponding to the paper "Solving a low-rank factorization model for matrix completion by a nonlinear successive over-relaxation algorithm", Math. Prog. Comp. (2012) 4: 333.