function visualize_errorcurves_combined(error_fro_rel,alg_name)
% visualize_errorcurves_combined(M,N,error_fro_rel,alg_name);

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
nr_algos=length(alg_name);
if ~isequal(nr_algos,length(error_fro_rel))
   error('Mismatch between provided number of algorithms and number of calculated error statistics.');
end

N=cell(1,nr_algos);
for i = 1:nr_algos
    N{i}=length(error_fro_rel{i});
end

maxIt=N{1};
minIt=N{1};
for l=1:nr_algos
   maxIt=max(maxIt,N{l});
   minIt=min(minIt,N{l});
end

figure
subplot(1,2,1)
for l=1:nr_algos
    plot(error_fro_rel{l},'Marker', 'o','MarkerSize',8);
    hold on
end
xlabel('iterations');
ylabel('Relative Frobenius error');
set(gca,'xtick',0:floor(maxIt./15):maxIt,'xticklabel',0:floor(maxIt./15):maxIt);
legend(alg_name,'FontSize',12);
hold on

subplot(1,2,2)
for l=1:nr_algos
   semilogy(error_fro_rel{l},'Marker', 'o','MarkerSize',8);
   hold on
end
xlabel('iterations');
ylabel('Logarithm of relative Frobenius error');
set(gca,'xtick',0:floor(maxIt./15):maxIt,'xticklabel',0:floor(maxIt./15):maxIt);
legend(alg_name,'FontSize',12);
title('Relative error vs. iteration count')
set(gcf,'numbertitle','off','name','Algorithmic comparisons for matrix completion') 
        
end

