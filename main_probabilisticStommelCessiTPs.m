clear; clc;
x=0:0.05:1.35;
y=0.09:0.05:1.69;
% for calculating nearest mesh point and assigning probabilities to mesh
% points
[X,Y] = meshgrid(x,y);
G = [X(:),Y(:)];
Frequency = zeros(length(G),1);
Frequency_FP = zeros(length(G),1);
% for strong AMOC state fold point
Frequency_FP2 = zeros(length(G),1);
etaSquared_cusp = 3;
%% etaSquared from UQLab results
load('UQpostSample_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES_steps400_NChains100.mat');
etaSquared = UQpostSample(1:1:end)';
etaSquared_synthData = 4;
M = length(etaSquared);
m_etaSquared = mean(etaSquared); % mean value of posterior sample

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

%% Define gray scale for plots
mygray = flipud(gray((length(Q)+5)/2)); % if length(Q) uneven otherwise adapt grayscale
mycolors = [mygray(2:end-1,:) ;flipud(mygray(2:end-2,:))];


%% Initialize vector of fold points that end the attracting part of the critical manifold
firstCritEtaP2 = find(etaSquared_sort>3,1);
numToDraw = N(length(Q))- max(firstCritEtaP2,N(1))+1;
x_crit = zeros(1,numToDraw);
y_crit = zeros(1,numToDraw);
% strong fold point
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
        x_crit(j-max(firstCritEtaP2,N(1))+1) = 2/3 +sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
        y_crit(j-max(firstCritEtaP2,N(1))+1) = h0(x_crit(j-max(firstCritEtaP2,N(1))+1),etaSquared_sort);
        
       %% strong fold point
        x_crit2(j-max(firstCritEtaP2,N(1))+1) = 2/3 -sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
        y_crit2(j-max(firstCritEtaP2,N(1))+1) = h0(x_crit2(j-max(firstCritEtaP2,N(1))+1),etaSquared_sort); 
        plot(y_crit(j-max(firstCritEtaP2,N(1))+1),x_crit(j-max(firstCritEtaP2,N(1))+1),'o','MarkerEdgeColor',mycolors(i,:),'MarkerFaceColor',mycolors(i,:),'MarkerSize',5);
    % plot tipping point strong AMOC state
        plot(y_crit2(j-max(firstCritEtaP2,N(1))+1),x_crit2(j-max(firstCritEtaP2,N(1))+1),'o','MarkerEdgeColor',mycolors(i,:),'MarkerFaceColor',mycolors(i,:),'MarkerSize',5);
    end
end
end

c=colorbar;
mapColors = mycolors(3:end,:);
colormap(mapColors);
set(c, 'YTick', linspace(1/6, 5/6, 3));
c.TickLabels = num2cell(massLevels);

xlabel('Critical nondim. freshwater flux');
ylabel('Critical salinity difference');
ax.FontSize = 12;
ax.Interpreter = 'latex';
xl=xlim;
xlim([0.8 xl(2)])
set(gca,'FontSize',12);

% % Save figure with probabilistic representation of tipping points in gray
% % scale
% savefig(fig1,strcat('probTPs_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.fig'));
% saveas(fig1,strcat('probTPs_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.eps'),'epsc');


%% Calculate fold points for etaSquare(1)=m_etaSquare
h0 = @(x,etaSquared) x.*(1+m_etaSquared*(1-x).^2);
y=h0(x,etaSquared);

%% Calculate fold points for etaSquared=val_syntheticData
h0 = @(x,etaSquared) x.*(1+etaSquared_synthData*(1-x).^2);
y=h0(x,etaSquared);

% https://stackoverflow.com/questions/54870941/finding-the-nearest-neighbor-to-a-single-point-in-matlab,
% last checked: 06.04.2021
[~,I] = pdist2(G, [x_crit' y_crit'], 'euclidean', 'Smallest', 1);
% from
% https://de.mathworks.com/matlabcentral/answers/142281-count-the-number-of-times-a-value-occurs-in-a-specific-of-an-array,
% last checked: 06.04.2021
edges = unique(I);
counts = histc(I, edges);
Frequency_FP(edges) = Frequency_FP(edges)+counts';

%% for strong AMOC state fold point
[~,I2] = pdist2(G, [x_crit2' y_crit2'], 'euclidean', 'Smallest', 1);
% from
% https://de.mathworks.com/matlabcentral/answers/142281-count-the-number-of-times-a-value-occurs-in-a-specific-of-an-array,
% last checked: 06.04.2021
edges2 = unique(I2);
counts2 = histc(I2, edges2);
Frequency_FP2(edges2) = Frequency_FP2(edges2)+counts2';

%% 2D histogram of location of critical values of nondimensional freshwater fluxes
fig2 = figure(2);
h1 = histogram(y_crit,'Normalization','pdf');
hold on
h2 = histogram(y_crit2,'Normalization','pdf');

leg = legend([h1 h2],{'$\mu^+$','$\mu^-$'});

leg.FontSize = 12;
leg.Interpreter = 'Latex';

xlabel('Critical nondim. freshwater flux');
ylabel('PDF estimate');
ax.FontSize = 12;
ax.Interpreter = 'latex';
xl=xlim;
xlim([0.85 xl(2)])
set(gca,'FontSize',12);

%% 2D histogram of location of critical values of salinity differences
fig3 = figure(3);
h1 = histogram(x_crit,'Normalization','pdf');
hold on
h2 = histogram(x_crit2,'Normalization','pdf');

leg = legend([h1 h2],{'$x^+$','$x^-$'});

leg.FontSize = 12;
leg.Interpreter = 'Latex';
leg.Location = 'north';

xlabel('Critical salinity difference');
ylabel('PDF estimate');
ax.FontSize = 12;
ax.Interpreter = 'latex';
xl=xlim;
xlim([0.4 xl(2)])
set(gca,'FontSize',12);

% %% Save 2D histograms on location of strong and weak TPs
% savefig(fig2,strcat('hist2DcritFreshwater_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.fig'));
% saveas(fig2,strcat('hist2DcritFreshwater_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.eps'),'epsc');
% 
% savefig(fig3,strcat('hist2DcritSalDiff_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.fig'));
% saveas(fig3,strcat('hist2DcritSalDiff_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.eps'),'epsc');