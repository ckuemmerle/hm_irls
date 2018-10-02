%% Script to reproduce experiment of Figure 5 of the paper:
% C. Kuemmerle, J. Sigl.
% "Harmonic Mean Iteratively Reweighted Least Squares for Low Rank Matrix
% Recovery", to appear in the Journal of Machine Learning Research, 2018.
% Available online: https://arxiv.org/abs/1703.05038v2.
rng(1e2*7)
%% Set parameters
% Number of rows of the matrix that defines the problem:
d1 = 40; 
% Number of columns:
d2 = 40; 
% Rank:
r = 10;
% Number of degrees of freedom of the setting:
df_LR = r*(d1 + d2 - r);
% Oversampling factor 
rho = 1.00;
% Number of measurements:
m = floor(min(rho*df_LR,d1*d2));
% % Maximal number of iterations 
opts_custom.N0 = 100;
% % (Non-convexity) parameter p (between 0 and 1)
opts_custom.p  = [0.8,0.5,0.25,0.01];
%% Sample matrix completion measurement matrix Phi
max_nr_resample = 5000;
[Phi,Omega] = sample_phi_MatrixCompletion(d1,d2,m,'resample',r,max_nr_resample);
% [aux_ind_Phi,nr_entr_row,nr_entr_col] = SamplePhi_MatrixCompletion_concise(d1,d2,m,r,nr_max_resample);
%% Sample the ground truth matrix X0
complexflag = 0;
modeX0      = 1 ;
[U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complexflag);
X0 = {U0,V0};
X0_full=U0*V0';
y = X0_full(Omega);

%% Choose algorithms for matrix completion
alg_name={'HM-IRLS','AM-IRLS','IRLS-col'};

%% Run algorithms for matrix completion
[Xr,outs,alg_name] = run_MC_algos(Phi,y,r,alg_name,opts_custom);
%% Calculating the error statistics
verbose = 1; % provide text output after error calculations
[error_fro_rel,error_fro] =get_frob_errors(Xr,X0,alg_name,verbose);
%% Visualization
visualize = 1;
if visualize == 1
    visualize_errorcurves_IRLScompare(error_fro_rel,alg_name,length(opts_custom.p));
end