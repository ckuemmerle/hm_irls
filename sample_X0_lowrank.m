function [U0,V0] = sample_X0_lowrank(d1,d2,r,modeX0,complex)
%sample_X0_lowrank This function constructs a random rank-r matrix
% X0 = U0 * V0' of dimension (d1 x d2). The randomness comes from Gaussian
% matrices of rank-r factor matrices.
%               Input:
%               d1      = number of rows of low-rank matrix X0 to be
%                           constructed
%           	d2      = number of columns of low-rank matrix X0  to be
%                           constructed
%               r       = rank of X0
%               modeX0  = 1: random model as in the paper [Kuemmerle, Sigl
%                           2018] (Section 5, X0 = U*S*V' where U,V are
%                           rank-r with i.i.d. standard Gaussian entries
%                           and S is (r x r), diagonal and has i.i.d.
%                           standard Gaussian entries.
%                       = 2: random model such that X0 = U*V' where U,V are
%                           rank-r with i.i.d. standard Gaussian entries.
%               complex = 1 for complex matrix X0 and 
%                         0 for real matrix X0.
%               Output: 
%               U0      = left factor of rank-r matrix X0
%               V0      = right factor of rank-r matrix X0
%
% The paper [Kuemmerle, Sigl 2018] can be found in:
% - "Harmonic Mean Iteratively Reweighted Least Squares for Low-Rank Matrix Recovery", 
%   to appear in the Journal of Machine Learning Research (JMLR).
%   Available online: https://arxiv.org/abs/1703.05038

if complex == 1
    U = (1/sqrt(2)).*(randn(d1,r)+1i.*randn(d1,r));
    V = (1/sqrt(2)).*(randn(d2,r)+1i.*randn(d2,r));
    if modeX0 == 1
        S = (1/sqrt(2)).*(randn(1,r)+1i.*randn(1,r));
    elseif modeX0 == 2
        S = ones(1,r);
    end
else
    U = randn(d1,r);
    V = randn(d2,r);
    if modeX0 == 1
        S = randn(1,r);
    elseif modeX0 == 2
        S = ones(1,r);
    end
end
U0 = U.*S;
V0 = V;

end

