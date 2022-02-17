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

% for nomenclature
m_etaSquaredval = strcat('m_etaP2',num2str(m_etaSquared));
m_etaSquaredval = strrep(m_etaSquaredval,'.','K');

%% prior sample of etaSquared
etaSquared_prior = (12.3-0.6)*rand(1,M) + 0.6;

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

%% Define gray scale for plots
mygray = flipud(gray((length(Q)+5)/2)); % if length(Q) uneven otherwise adapt grayscale
mycolors = [mygray(2:end-1,:) ;flipud(mygray(2:end-2,:))];


%% Initialize vector of fold points that limit the attracting part of the critical manifold
firstCritEtaP2 = find(etaSquared_sort>3,1);
% %% comment if if firstCritEtaP2>N(1)
% numToDraw = N(length(Q))- N(1);
% uncomment if if firstCritEtaP2>N(1)
numToDraw = N(length(Q))- max(firstCritEtaP2,N(1))+1;
x_crit = zeros(1,numToDraw);
y_crit = zeros(1,numToDraw);
% strong fold point
x_crit2 = zeros(1,numToDraw);
y_crit2 = zeros(1,numToDraw);
distBif = cell(1,numToDraw);

%% Calculate fold points that end attracting part of the critical manifold and plot critical manifolds
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
    % plot tipping component freshwater flux weak AMOC state
        plot(etaSquared_sort(j),y_crit(j-max(firstCritEtaP2,N(1))+1),'o','MarkerEdgeColor',mycolors(i,:),'MarkerFaceColor',mycolors(i,:),'MarkerSize',5);
    % plot tipping component freshwater flux strong AMOC state
        plot(etaSquared_sort(j),y_crit2(j-max(firstCritEtaP2,N(1))+1),'o','MarkerEdgeColor',mycolors(i,:),'MarkerFaceColor',mycolors(i,:),'MarkerSize',5);
    end
end
end

plot([3 3], ylim, 'r', 'LineWidth', 2);

c=colorbar;
mapColors = mycolors(3:end,:);
colormap(mapColors);
set(c, 'YTick', linspace(1/6, 5/6, 3));
c.TickLabels = num2cell(massLevels);

xlabel('\eta^2');
ylabel('Critical nondim. freshwater flux');
ax.FontSize = 12;
ax.Interpreter = 'latex';
xl=xlim;
xlim([2.5 xl(2)])
set(gca,'FontSize',12);

% %% Save figure
% savefig(fig1,strcat('probCritFreshwater_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.fig'));
% saveas(fig1,strcat('probCritFreshwater_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.eps'),'epsc');