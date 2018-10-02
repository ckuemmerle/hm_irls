% Script to illustrate the example of Figure 2 of the paper:
% "Harmonic Mean Iteratively Reweighted Least Squares for Low Rank Matrix
% Recovery", to appear in the Journal of Machine Learning Research, 2018.
% Available online: https://arxiv.org/abs/1703.05038v2:
% Comparison of different IRLS variants for low-rank matrix competlion.
%% Dimensions of matrix X0 to be recovered:
d1= 4;  % number of rows
d2= 4;  % number of columns
r = 1;  % rank
%% Definition of mask which indicates the known indices
Phi=[0,0,0,1;1,0,0,1;0,1,0,0;1,1,1,0];
Omega= find(Phi);
m = length(Phi);
%% Definition of matrix X0
U0=[1e0,1e1,-2e0,1e-1]'; Sigma = 1;
% %  U0=U0./norm(U0);
V0=[1,2,3,4]';
X0      = {U0*Sigma,V0};
X0_full = U0*Sigma*V0';
y       = X0_full(Omega);

U0_norm=U0./norm(U0);
P_U=U0_norm*U0_norm';
V0_norm=V0./norm(V0);
P_V=V0_norm*V0_norm';
coh_U=max(d1./r.*[norm(P_U*[1;0;0;0])^2;norm(P_U*[0;1;0;0])^2;norm(P_U*[0;0;1;0])^2;norm(P_U*[0;0;0;1])^2]);
coh_V=max(d2./r.*[norm(P_V*[1;0;0;0])^2;norm(P_V*[0;1;0;0])^2;norm(P_V*[0;0;1;0])^2;norm(P_V*[0;0;0;1])^2]);

%% Choose algorithms for matrix completion
alg_name={'HM-IRLS','IRLS-col','IRLS-row','AM-IRLS'};
opts_custom.p=[0.1]; %%% choose (non-)convexity parameters p for IRLS

%% Run algorithms for matrix completion
[Xr,outs,alg_name] = run_MC_algos(Phi,y,r,alg_name,opts_custom);
nr_algos=length(alg_name);
%% Calculate the H matrices
H=cell(1,nr_algos);
for l=1:nr_algos
   singval_c=outs{l}.singval;
   for k=1:outs{l}.N
       H{l}{k}=zeros(d1,d2);
       for ii=1:d1
           for jj=1:d2
               if l == 1
                    H{l}{k}(ii,jj)=2.*((singval_c{k}(ii)^2+outs{l}.eps(k).^2)^((2-opts_custom.p)/2)+(singval_c{k}(jj)^2+outs{l}.eps(k).^2)^((2-opts_custom.p)/2))^(-1);
               elseif l == 2
                    H{l}{k}(ii,jj)=(singval_c{k}(ii)^2+outs{l}.eps(k).^2)^((opts_custom.p-2)/2);
               elseif l == 3
                    H{l}{k}(ii,jj)=(singval_c{k}(jj)^2+outs{l}.eps(k).^2)^((opts_custom.p-2)/2);
               elseif l == 4
                    H{l}{k}(ii,jj)=0.5.*((singval_c{k}(ii)^2+outs{l}.eps(k).^2)^((opts_custom.p-2)/2)+(singval_c{k}(jj)^2+outs{l}.eps(k).^2)^((opts_custom.p-2)/2));
               end
           end
       end
   end
end

%% Create figure of H matrices after first iteration
figure
for i=1:nr_algos
    subplot(1,nr_algos,i)
    imagesc(H{i}{1});
    set(gca,'ytick',[1 2 3 4]);
    if i == 1
        xlabel('Column j');
        ylabel('Row i');
    end
end
 colbar = colorbar;
% colorbar('off');
%% Calculating the error statistics
verbose = 1; % provide text output after error calculations
[error_fro_rel,error_fro] =get_frob_errors(Xr,X0,alg_name,verbose);
%% Visualization of Frobenius errors of iterates of algorithms
visualize = 1;
if visualize == 1
    visualize_errorcurves_combined(error_fro_rel,alg_name);
end