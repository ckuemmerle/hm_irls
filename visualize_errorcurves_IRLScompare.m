function visualize_errorcurves_IRLScompare(error_f_rel,alg_name,sizep)
%visualize_errorcurves_IRLScompare Plots a semilogarithmic plot of
%(relative) Frobenius errors given in error_f_rel.
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
nr_algos=length(alg_name);
maxIt=length(error_f_rel{1});
minIt=length(error_f_rel{1});%N{1};
colorscheme= hsv(sizep);
cmap=zeros(sizep,3);
cmap(4,:)=colorscheme(sizep,:);
cmap(3,:)=[0 0 0];
cmap(2,:)=[1 0 0];
cmap(1,:)=[0 1 0];
for l=1:nr_algos
   maxIt=max(maxIt,length(error_f_rel{l}));
   minIt=min(minIt,length(error_f_rel{l}));
end
alg_name_disp=cell(1,nr_algos);

figure
ii = 0;
for j=[0,1:1:sizep-1]
    jj = j;
    if j == 0
        jj = sizep;
    end
    for l=1:nr_algos
       if mod(l-1,sizep) == j
           ii=ii+1;
           if floor((l-1)./sizep) == 0
               linetype = '-';
               markertype = 'o';
               semilogy(error_f_rel{l},'Marker', markertype,'LineStyle',linetype,'Color',cmap(jj,:),'MarkerSize',5,'LineWidth',2);
               ylim([1e-10,1e0])
               hold on
           elseif floor((l-1)./sizep) == 1
               linetype = '--';
               semilogy(error_f_rel{l},'LineStyle',linetype,'Color',cmap(jj,:),'LineWidth',2);
               ylim([1e-10,1e0])
               hold on
           elseif floor((l-1)./sizep) == 2
               linetype = ':';
               semilogy(error_f_rel{l},'LineStyle',linetype,'Color',cmap(jj,:),'LineWidth',2);
%                ylim([1e-2,1e0])
               hold on
           end
            alg_name_disp{ii}=alg_name{l};
       end
    end
end

xlabel('Iteration n','FontSize',15);
ylabel('Relative Frobenius error $\|X^{(n)}-X_0\|_F/\|X_0\|_F$','Interpreter','latex','FontSize',15);
yticks([10^(-2),10^(-1),10^(0)])
ax=gca;
ax.FontSize=14;
legend(alg_name_disp,'FontSize',17);
hold off
set(gcf,'numbertitle','off','name','Algorithmic comparisons for matrix completion') 

        
end

