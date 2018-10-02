function [error_fro_rel,error_fro] =get_frob_errors(Xr,X0,alg_name,verbose)

nr_algos=length(Xr);
error_fro       = cell(1,nr_algos);
error_fro_rel   = cell(1,nr_algos);
for l=1:nr_algos
    error_fro_rel{l}= zeros(1,length(Xr{l}));
    error_fro{l}    = zeros(1,length(Xr{l}));
    for i = 1:length(Xr{l})
        [error_fro_rel{l}(i),error_fro{l}(i)] = ...
        calc_frob_error(Xr{l}{i},X0);
    end
    if verbose
        disp(alg_name{l});
        disp(['Relative Frob. error to X0:      ',num2str(error_fro_rel{l}(length(Xr{l})))]);
    end
    
    
end
end

function [error_fro_rel,error_fro] = ...
    calc_frob_error(Xr_c,X0)

if iscell(Xr_c)
    if size(Xr_c,2) == 2
        U=Xr_c{1};
        V=Xr_c{2};
        
        error_fro = trace((U'*U)*(V'*V));
        if iscell(X0)
            U0=X0{1};
            V0=X0{2};
            error_fro = - 2* real(trace((U0'*U)*(V'*V0)))         + error_fro;
            error_fro = + trace((U0'*U0)*(V0'*V0))          + error_fro;
        else
            error_fro = - 2* real(trace(V'*X0*U))                 + error_fro;
            error_fro = + norm(X0,'fro').^2                 + error_fro;
        end
        error_fro =   sqrt(error_fro);
    else
        error('Error in the format of the recovered matrices.')
        
    end
else
    if iscell(X0)
        U0=X0{1};
        V0=X0{2};
        X0 = U0*V0';
    end
    error_fro = norm(Xr_c-X0,'fro');
end
error_fro=real(error_fro);
if iscell(X0)
    error_fro_rel = error_fro./sqrt(real(trace((U0'*U0)*(V0'*V0))));
else
    error_fro_rel = error_fro./norm(X0,'fro');
end

end