function opts = getDefaultOpts_IRLS

opts.p          = 0.01;  % Schatten-p parameter (p close to 0 for very non-convex objective, p=1 for convex)
opts.N0         = 2000;  % max. number of iteration
opts.tol        = 1e-7;  % stopping criterion, lower bound on relative change of Frobenius norm
opts.epsmin     = 1e-9;  % minimal value for epsilon smoothing
opts.linsolve   = 'direct_solve';

end

