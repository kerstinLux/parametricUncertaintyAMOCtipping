clear; clc;
%% Attention: graphic intense calculation, otherwise reduce number of samples M from posterior used
%% Plot critical manifold for sample of values of the ratio etaSquared of the diffusive time scale to the advective one (Cessi(1994))
x=0:0.001:1.35;
y=0.09:0.001:1.69;
[X,Y] = meshgrid(x,y);
G = [X(:),Y(:)];
Frequency = zeros(length(G),1);
Frequency_FP = zeros(length(G),1);
% for strong AMOC state fold point
Frequency_FP2 = zeros(length(G),1);
etaSquared_cusp = 3; % below, there is no fold bifurcation
%% load etaSquared from UQLab results
load('UQpostSample_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES_steps400_NChains100.mat');
etaSquared = UQpostSample(1:10:end)'; % use full posterior sample loaded
m_etaSquared = mean(etaSquared);
M = length(etaSquared);
% %% sample etaSquared from random normal sample
% m_etaSquared = 4;
% std_etaSquared = 1;
% M=10^3;
% etaSquared = m_etaSquared + std_etaSquared*randn(1,M);

m_etaSquaredval = strcat('m_etaP2',num2str(m_etaSquared));
m_etaSquaredval = strrep(m_etaSquaredval,'.','K');

%% Define mass levels for gray color code
massLevels = [0.50 0.75 0.975]; %ATTENTION: if adaptation desired, remember to also adapt corresponding quantile levels "levels"
%% Note that the quantile levels are constructed such that the dark gray area comprises 50% probability mass,
%  the middle gray area inclduing the dark gray one comprises 75% and the
% whole gray area 97.5% --> corresponds to the standard definition of a
% confidence interval
% only the symmetric outer 2.5% of etaSquared realizations are not plotted/shown
levels = [0.0125 0.125 0.25 0.75 0.875 0.9875]; % quantile levels needed for confLevels = [0.025 0.25 0.50 0.75 0.975];
Q = quantile(etaSquared,levels);
etaSquared_sort = sort(etaSquared);
N = zeros(size(Q));
for i=1:length(Q)
[~,N(i)]=min(abs(etaSquared_sort-Q(i)));
end
%% Determine number of realizations with fold points
firstCritEtaP2 = find(etaSquared_sort>3,1);
numToDraw = N(length(Q))- firstCritEtaP2+1; % #realizations with fold points

%% Define gray scale for plots
mygray = flipud(gray((length(Q)+5)/2)); % if length(Q) uneven otherwise adapt grayscale
mycolors = [mygray(2:end-1,:) ;flipud(mygray(2:end-2,:))];


%% Initialize vector of fold points that end the attracting part of the critical manifold
% for weak fold point
x_crit = zeros(1,numToDraw);
y_crit = zeros(1,numToDraw);
% for strong fold point
x_crit2 = zeros(1,numToDraw);
y_crit2 = zeros(1,numToDraw);
distBif = cell(1,numToDraw);

%% Calculate fold points that limit attracting part of the critical manifold and plot critical manifolds
fig1 = figure(1);
hold on
for i=1:length(Q)-1
for j = N(i)+1:N(i+1)
    h0 = @(x,etaSquared) x.*(1+etaSquared(j)*(1-x).^2);
    y=h0(x,etaSquared_sort);
    l1 = plot(y,x,'color',mycolors(i,:),'LineWidth',2);
    hold on
%     %% comment if firstCritEtaP2>N(1)
%         %% Critical manifold splits into three parts: calculate fold point where attracting part ends
%         x_crit(j-N(1)) = 2/3 +sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
%         y_crit(j-N(1)) = h0(x_crit(j-N(1)),etaSquared_sort);
%         
%        %% strong fold point
%        x_crit2(j-N(1)) = 2/3 -sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
%         y_crit2(j-N(1)) = h0(x_crit2(j-N(1)),etaSquared_sort);   
% uncomment if if firstCritEtaP2>N(1)
    if etaSquared_sort(j)>3
        %% Critical manifold splits into three parts: calculate fold point where attracting part ends
        x_crit(j-firstCritEtaP2+1) = 2/3 +sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
        y_crit(j-firstCritEtaP2+1) = h0(x_crit(j-firstCritEtaP2+1),etaSquared_sort);
        
       %% strong fold point
       x_crit2(j-firstCritEtaP2+1) = 2/3 -sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
        y_crit2(j-firstCritEtaP2+1) = h0(x_crit2(j-firstCritEtaP2+1),etaSquared_sort);
    end
end
end

%% Calculate fold points for etaSquare(1)=m_etaSquare
h0 = @(x,etaSquared) x.*(1+m_etaSquared*(1-x).^2);
y=h0(x,etaSquared);
l1 = plot(y,x,'k','LineWidth',2);
hold on
%% Calculate critical manifold for etaSquare(1)=etaSquared_cusp
h0 = @(x,etaSquared) x.*(1+etaSquared_cusp*(1-x).^2);
y=h0(x,etaSquared);
l1 = plot(y,x,'r','LineWidth',2);
hold on

% inspired from https://stackoverflow.com/questions/54870941/finding-the-nearest-neighbor-to-a-single-point-in-matlab,
% last checked: 06.04.2021
%% for weak AMOC state fold point
[~,I] = pdist2(G, [x_crit' y_crit'], 'euclidean', 'Smallest', 1);
% inspired from https://de.mathworks.com/matlabcentral/answers/142281-count-the-number-of-times-a-value-occurs-in-a-specific-of-an-array,
% last checked: 06.04.2021
edges = unique(I);
counts = histc(I, edges);
Frequency_FP(edges) = Frequency_FP(edges)+counts';
%% for strong AMOC state fold point
[~,I2] = pdist2(G, [x_crit2' y_crit2'], 'euclidean', 'Smallest', 1);
% inspired from https://de.mathworks.com/matlabcentral/answers/142281-count-the-number-of-times-a-value-occurs-in-a-specific-of-an-array,
% last checked: 06.04.2021
edges2 = unique(I2);
counts2 = histc(I2, edges2);
Frequency_FP2(edges2) = Frequency_FP2(edges2)+counts2';

xlabel('Nondimensional freshwater flux');
ylabel('Salinity difference');
ax.FontSize = 12;
ax.Interpreter = 'latex';
set(gca,'FontSize',12);

c=colorbar;
mapColors = mycolors(3:end,:);
colormap(mapColors);
set(c, 'YTick', linspace(0.125, 0.875, 3));
c.TickLabels = num2cell(massLevels);


% %% Save figure
% % for etaSquared from posterior distribution resulting from UQLab Bayesian
% % inversion
% savefig(fig1,strcat('probBifCurves_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3','.fig'));
% saveas(fig1,strcat('probBifCurves_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3','.eps'),'epsc');
% % %% for etaSquared from posterior distribution resulting from random normal sample
% % savefig(fig1,strcat('completeProbBifDiagram_Stommel_normal_mean_4_sigma1','.fig'));
% % saveas(fig1,strcat('completeProbBifDiagram_Stommel_normal_mean_4_sigma1','.eps'),'epsc');