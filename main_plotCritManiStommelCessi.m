clear; clc;
%% Plot critical manifold for various values of the ratio etaSquared of the diffusive time scale to the advective one (Cessi(1994))
fig1 = figure(1)
etaSquared = [1 3 4 7.5];
myColors = {'b','r','m','k'};
for i=1:length(etaSquared)-2
    h0 = @(x) x.*(1+etaSquared(i)*(1-x).^2);
    x = 0:0.01:1.4;
    y = h0(x);
    l(i) = plot(y,x,'Color',myColors{i},'LineWidth',2);
    hold on
end

for i=length(etaSquared)-1:length(etaSquared)
    %% weak AMOC fold point
    h0 = @(x) x.*(1+etaSquared(i)*(1-x).^2);
    x_crit = 2/3 +sqrt(4/9 - (1+etaSquared(i))/(3*etaSquared(i)));
    y_crit = h0(x_crit);     
    %% strong AMOC fold point
    x_crit2 = 2/3 -sqrt(4/9 - (1+etaSquared(i))/(3*etaSquared(i)));
    y_crit2 = h0(x_crit2);   
    x_att1 = x_crit:0.001:1.4;
    x_rep = x_crit2:0.001:x_crit;
    x_att2 = 0:0.001:x_crit2;
    y_att1=h0(x_att1);
    y_rep = h0(x_rep);
    y_att2 = h0(x_att2);
    l(i) = plot(y_att1,x_att1,'Color',myColors{i},'LineWidth',2);
    hold on
    plot(y_rep,x_rep,'Color',myColors{i},'LineStyle','--','LineWidth',2);
    plot(y_att2,x_att2,'Color',myColors{i},'LineWidth',2);

    plot(y_crit,x_crit,'o','Markersize',5,'MarkerFaceColor',myColors{i},'MarkerEdgeColor',myColors{i});
    plot(y_crit2,x_crit2,'o','Markersize',5,'MarkerFaceColor',myColors{i},'MarkerEdgeColor',myColors{i});
end

xlabel('Nondim. freshwater flux');
ylabel('Salinity difference');
ax = gca;
ax.FontSize = 12;
set(gca,'FontSize',12);
leg = legend([l(1) l(2) l(3) l(4)],'$\eta^2 = 1$','$\eta^2 = 3$','$\eta^2 = 4$','$\eta^2 = 7.5$');
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Southeast';

% savefig(fig1,strcat('redStommelCessi_critManifold.fig'));
% saveas(fig1,strcat('redStommelCessi_critManifold.eps'),'epsc');
