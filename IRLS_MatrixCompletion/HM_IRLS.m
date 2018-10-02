function [X,N,eps,singval,time] = ...
    HM_IRLS(prob,opts)
%HM-IRLS Harmonic mean iteratively least squares algorithm, implemented for
%matrix completion measurements.
%
%   Inputs: prob = Struct with fields "d1", "d2", "r", "Omega", "Phi" for
%                       problem speficiation
%           opts = Struct with fields "N0", "p", "epsmin", "tol" and
%                       "linsolve"
%  Outputs:    X = (1 x N) cell with (d1 x d2) matrices: 
%                       iterates of HM-IRLS algorithm
%              N = number of performed iterations
%            eps = (1 x N) vector with values of epsilon-smoothing at every
%                   iteration
%        singval = (1 x N) cell with vectors of singular values of X{l}
%           time = elapsed time
%
% For more information, see the paper:
% C. Kuemmerle, J. Sigl.
% "Harmonic Mean Iteratively Reweighted Least Squares for Low Rank Matrix
% Recovery", to appear in the Journal of Machine Learning Research, 2018.
% Available online: https://arxiv.org/abs/1703.05038.
%
% Written by Christian Kuemmerle, Technical University of Munich.
% Questions and bug reports to: c.kuemmerle@tum.de.
% Version 1.1, Updated 10/01/2018.

%% Get problem specification and algorithmic options
d1          = prob.d1;       % first dimension of matrix to be recovered
d2          = prob.d2;       % second dimension of matrix to be recovered
r           = prob.r;        % target rank
Omega       = prob.Omega;    % indicates "seen" linear indices of matrix
Phi         = prob.Phi;      % measurement mask
y           = sparse(prob.y);% measured entries

N0          = opts.N0;      % max. number of iterations
p           = opts.p;       % Schatten-p parameter
epsmin      = opts.epsmin;  % minimal eps
tol         = opts.tol;     % stopping criterion (minimal relative change of Frobenius norm)

dd=min(d1,d2);

%%% get some auxiliary quantities from the sampling mask
nr_entr_col = sum(Phi,1)';
nr_entr_row = sum(Phi,2);
[rowind_mask,colind_mask]   = find(Phi);
[~,ind_reordered] = sort(rowind_mask);
%% Iterations of the algorithm
X                = cell(1,N0);
eps              = zeros(1,N0);
singval          = cell(1,N0);

l=1;
W_2_1_inv=eye(d1);
W_2_2_inv=eye(d2);
X_c      =zeros(d1,d2);
eps_c    =Inf;
tic
while l <= N0
    X_cm1=X_c;
    %% Calculate the solution of the minimization problem
    if mod(l,25) == 0 
        disp(['Begin Iteration ',num2str(l),'...']);
    end
    
    if strcmp(opts.linsolve,'direct_solve')
        %% Update matrix of linear system to be solved
        PhWiPh_star = update_PhWiPh(W_2_1_inv,W_2_2_inv,rowind_mask,colind_mask,...
            ind_reordered,nr_entr_col,nr_entr_row);
        %% Calculate new iterate
            X_c = directsolve_weightedLS(PhWiPh_star,y,...
               W_2_1_inv,W_2_2_inv,Omega,'harmonic');
    end
    %% Calculate weight matrices W_L^{-1} and W_R^{-1} from SVD information
    [U_X_c,singval_mat_c,V_X_c] = svd(X_c,'econ');
    singval_c = diag(singval_mat_c);
    eps_c = min(eps_c,max(singval_c(r+1),epsmin));

    Sigma_c_diag= (singval_c.^2+ eps_c.^2.*ones(dd,1)).^((2-p)/2);
    U_X_c_Sig= U_X_c*spdiags(Sigma_c_diag,0,speye(dd,dd));
    V_X_c_Sig= V_X_c*spdiags(Sigma_c_diag,0,speye(dd,dd));
    W_2_1_inv = U_X_c_Sig*U_X_c';
    W_2_2_inv = V_X_c_Sig*V_X_c'; 
    
    eps(l)      = eps_c;
    singval{l}  = singval_c;
    X{l}        = X_c;
    %% Check stopping criterion
    if l >= 2 && l <= N0 && norm(X_c-X_cm1,'fro')/norm(X_c,'fro') < tol
        N0 = l;
    end
    l=l+1;
end
time=toc;
%% Prepare output
X             = X(1:N0);
eps           = eps(1,1:N0);
singval       = singval(1,1:N0);
N             = N0;
end





function PhWiPh_star = update_PhWiPh(W_2_1_inv,W_2_2_inv,rowind_mask,colind_mask,...
            ind_reordered,nr_entr_col,nr_entr_row)
d1=size(W_2_1_inv,1);
d2=size(W_2_2_inv,1);
m =length(rowind_mask);
%% Calculate sparse matrix Phi*(Id,kron(W_L^(-1)))*(Phi^*) 
ind=1; count=0;
ind_row=zeros(floor(m*m./d2)+1,1);
ind_col=zeros(floor(m*m./d2)+1,1);
vec_val=zeros(floor(m*m./d2)+1,1);
for i=1:d2
    val = W_2_1_inv(rowind_mask(ind-1+(1:nr_entr_col(i))),rowind_mask(ind-1+(1:nr_entr_col(i))));
    [ind_i,ind_j] = ind2sub([nr_entr_col(i),nr_entr_col(i)],reshape(1:(nr_entr_col(i)^2),[nr_entr_col(i),nr_entr_col(i)]));
    ind_row(count+(1:(nr_entr_col(i)^2))) =  ind-1+ind_i;
    ind_col(count+(1:(nr_entr_col(i)^2))) =  ind-1+ind_j;
    vec_val(count+(1:(nr_entr_col(i)^2))) =  val;

    count = count+nr_entr_col(i)^2;
    ind=ind+nr_entr_col(i);
end
ind_row=ind_row(1:count);
ind_col=ind_col(1:count);
vec_val=vec_val(1:count);
PhWiPh_star_l = sparse(ind_row,ind_col,vec_val,m,m); % calculate Phi*(kron(Id,W_L^(-1)))*(Phi^*) 
%% Calculate sparse matrix Phi*(kron(W_R^(-1),Id))*(Phi^*)
ind=1; count=0;
ind_row=zeros(floor(m*m./d1)+1,1);
ind_col=zeros(floor(m*m./d1)+1,1);
vec_val=zeros(floor(m*m./d1)+1,1);
for i=1:d1
    val = W_2_2_inv(colind_mask(ind_reordered(ind-1+(1:nr_entr_row(i)))),colind_mask(ind_reordered(ind-1+(1:nr_entr_row(i)))));
    [ind_i,ind_j] = ind2sub([nr_entr_row(i),nr_entr_row(i)],reshape(1:(nr_entr_row(i)^2),[nr_entr_row(i),nr_entr_row(i)]));
    ind_row(count+(1:(nr_entr_row(i)^2))) =  ind_reordered(ind-1+ind_i);
    ind_col(count+(1:(nr_entr_row(i)^2))) =  ind_reordered(ind-1+ind_j);
    vec_val(count+(1:(nr_entr_row(i)^2))) =  val;
    count = count+nr_entr_row(i)^2;
    ind=ind+nr_entr_row(i);
end
ind_row=ind_row(1:count);
ind_col=ind_col(1:count);
vec_val=vec_val(1:count);
PhWiPh_star_r = sparse(ind_row,ind_col,vec_val,m,m); % calculate Phi*(kron(W_R^(-1),Id))*(Phi^*) 

PhWiPh_star = 0.5.*(PhWiPh_star_l+PhWiPh_star_r); % calculate Phi*\tilde(W)^{-1}*(Phi^*) 
end