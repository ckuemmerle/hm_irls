% Script to compare the algorithm of the paper [1] with 
% - "LMaFit" (Zaiwen Wen, Wotao Yin and Yin Zhang 2012), see [1] for full
% reference,
% - "LRGeomCG" (Bart Vandereycken 2013), see [1] for full reference,
% on random matrix completion examples.
%
% [1] "Harmonic Mean Iteratively Reweighted Least Squares for Low Rank Matrix
% Recovery", to appear in the Journal of Machine Learning Research, 2018.
% Available online: https://arxiv.org/abs/1703.05038v2.
%% Set parameters
% Number of rows of the matrix that defines the problem:
d1 = 100; 
% Number of columns:
d2 = 100; 
% Rank:
r = 7;
% Number of degrees of freedom of the setting:
df_LR = r*(d1 + d2 - r);
% Number of measurements:
m = floor(min(1.05*df_LR,d1*d2));
%% Sample the measurement matrix Phi (pattern of revealed entries) 
max_nr_resample = 1000;
[Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
 
%% Sample the ground truth matrix X0 and measured entries y
modeX0      = 1;
complexflag = 0;
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag);
X0 = {U0,V0};
X0_full=U0*V0';
y=X0_full(Omega);

%% Choose algorithms for matrix completion
alg_name={'HM-IRLS','LMaFit','LRGeomCG'};
opts_custom.p=[0.1,0.5]; %%% choose (non-)convexity parameters p for IRLS

%% Run algorithms for matrix completion
[Xr,outs,alg_name] = run_MC_algos(Phi,y,r,alg_name,opts_custom);

%% Calculating the error statistics
verbose = 1; % provide text output after error calculations
[error_fro_rel,error_fro] =get_frob_errors(Xr,X0,alg_name,verbose);
%% Visualization of Frobenius errors of iterates of algorithms
visualize = 1;
if visualize == 1
    visualize_errorcurves_combined(error_fro_rel,alg_name);
end