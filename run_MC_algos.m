function [Xr,outs,alg_name] = run_MC_algos(Phi,y,r,alg_name,opts_custom)
%run_MC_algos This function runs different algorithms (indicated by alg_name)
%for a given matrix completion problem with entry mask Phi and provided data y.

%       Inputs: 
%             Phi = (d1 x d2) sparse matrix. Indicating the mask
%                       corresponding to known indices.
%               y = (m x 1) vector of samples. Use NaN for sample entries
%                       that are not known (completion case).
%               r = target rank to be used
%        alg_name = (1 x nr_algos) cell of character strings indicating
%                    which algorithm to use.
%     opts_custom = (1 x nr_algos) cell of custom options for the
%                   respective algorithms.


[d1,d2]=size(Phi);
Omega = find(Phi);
nr_algos = length(alg_name);
opts     = cell(1,nr_algos);
opts_new = cell(1,nr_algos);

%%% prepare multiple algorithm instances of IRLS variants if multiple
%%% values for quasi-norm parameter p are provided
if isfield(opts_custom,'p')
    len_p    = length(opts_custom.p);
    if len_p >= 1
        ind=find(contains(alg_name,'HM-IRLS')+contains(alg_name,'AM-IRLS')+...
            contains(alg_name,'IRLS-col')+contains(alg_name,'IRLS-row'));
        for i=1:length(ind)
            ind_c=ind(i)+(i-1)*(len_p-1);
            alg_name_ind_c=cell(1,len_p);
            for ii=1:len_p
                alg_name_ind_c{ii}=[alg_name{ind_c},' p=',num2str(opts_custom.p(ii))];
            end
            if ind_c < length(alg_name)
                alg_name= [alg_name(1:ind_c-1),alg_name_ind_c,alg_name(ind_c+1:end)];
                opts_new= [opts_new(1:ind_c-1),repelem(opts_new(ind_c),len_p),opts_new(ind_c+1:end)];
            else
                alg_name= [alg_name(1:ind_c-1),alg_name_ind_c];
                opts_new= [opts_new(1:ind_c-1),repelem(opts_new(ind_c),len_p)];
            end
            for j=1:len_p
                opts_new{ind_c+j-1}.p=opts_custom.p(j);
            end
        end
    end
end
nr_algos = length(alg_name);
names = fieldnames(opts_custom);
for l=1:nr_algos
    for i=1:length(names)
        if not(isequal(names(i),{'p'}))
            opts_new{l}.(names{i})=opts_custom.(names{i});
        end
    end
end



prob  = cell(1,nr_algos);
outs  = cell(1,nr_algos);
Xr    = cell(1,nr_algos);

%%% run algorithm (with default options, priority on option in opts_new if
%%% applicable)
for l=1:nr_algos
    if contains(alg_name{l},'HM-IRLS') || contains(alg_name{l},'AM-IRLS') ...
            || contains(alg_name{l},'IRLS-col') || contains(alg_name{l},'IRLS-row')
        opts{l} = getDefaultOpts_IRLS(alg_name{l});
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        prob{l}.d1     = d1;
        prob{l}.d2     = d2;
        prob{l}.r      = r;
        prob{l}.Omega  = Omega;
        prob{l}.Phi    = Phi;
        prob{l}.y      = y;
        if contains(alg_name{l},'HM-IRLS')
                [Xr{l},N,eps,singval,time] = ...
            HM_IRLS(prob{l},opts{l});
        elseif contains(alg_name{l},'AM-IRLS')
                [Xr{l},N,eps,singval,time] = ...
            AM_IRLS(prob{l},opts{l});
        elseif contains(alg_name{l},'IRLS-col')
            opts{l}.variant='left';
                [Xr{l},N,eps,singval,time] = ...
            IRLS_onesided(prob{l},opts{l});
        elseif contains(alg_name{l},'IRLS-row')
            opts{l}.variant='right';
                [Xr{l},N,eps,singval,time] = ...
            IRLS_onesided(prob{l},opts{l});  
        end
        outs{l}.N   = N;
        outs{l}.eps = eps;
        outs{l}.singval = singval;
        outs{l}.time    = time;
        
    elseif strcmp(alg_name{l},'LMaFit')
        opts{l}.tol=1.25e-10;
        opts{l}.print=0;
        opts{l}.N0=500; 
        opts{l}.est_rank=0; 
        opts{l}.Zfull=1;
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l}.maxit=opts{l}.N0;
        tic
        [~,~,Out] = lmafit_mc_adp(d1,d2,r,Omega,y,opts{l});
        outs{l}.N=Out.iter-1;
        outs{l}.time=toc;
        for kk=1:outs{l}.N
            Xr{l}{kk}={Out.Xhist{kk},Out.Yhist{kk}'};
        end
       
   elseif strcmp(alg_name{l},'LRGeomCG')
        tic
        prob{l} = make_prob_Riem(y,Omega,d1,d2,r);
        start = make_start_x_Riem(prob{l});
        opts{l}.N0 = 500;
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        opts{l} = default_opts_Riem(opts{l}.N0);
        opts{l} = setExtraOpts(opts{l},opts_new{l});
        [~,~,~,Xr{l},outs{l}.N] = LRGeomCG_outp(prob{l},opts{l},start);
        outs{l}.time=toc;

%    elseif strcmp(alg_name{l},'NNM')
%        tic
% %        [~,outs{l}.x_evol] = hankel_NNM_cvx(y(Omega).',Omega.',floor(n/2),n);
%         yr{l}=outs{l}.x_evol.';
%         outs{l}.time=toc;
   end
end


end

function opts = setExtraOpts(opts,opts_new)
    if ~isempty(opts_new)
        optNames = fieldnames(opts_new);
        for i = 1:length(optNames)
            optName = optNames{i};
            opts.(optName) = opts_new.(optName);
        end
    end

end

