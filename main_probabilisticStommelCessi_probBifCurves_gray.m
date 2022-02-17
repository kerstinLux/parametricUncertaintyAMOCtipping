clear; clc;
%% Plot critical manifold for various values of the ratio etaSquared of the diffusive time scale to the advective one (Cessi(1994))
x=0:0.001:1.35;
% for calculating nearest mesh point and assigning probabilities to mesh
% points
y=0.09:0.001:1.69;
[X,Y] = meshgrid(x,y);
G = [X(:),Y(:)];
Frequency = zeros(length(G),1);
Frequency_FP = zeros(length(G),1);
% for strong AMOC state fold point
Frequency_FP2 = zeros(length(G),1);
etaSquared_cusp = 3;

%% etaSquared from UQLab results % comment if choice of standard distribution desired
load('UQpostSample_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES_steps400_NChains100.mat');
etaSquared = UQpostSample(1:10:end)'; % use every 10th posterior sample loaded
M = length(etaSquared); % number of posterior samples
m_etaSquared = mean(etaSquared); % mean value of posterior sample used from UQLab results
etaSquared_synthData = 4; % to adapt according to loaded data!



%% Uncomment if a given distribution shall be used and comment load command
%% etaSquared from choice of standard distribution
% %% bounds for etaSquared, the value of ratio of diffusive and advective timescale
% a = 0.6; b = 12.3;
% %% Simulate M etaSquare realizations of etaSquared
% M = 10^3;
% etaSquared = zeros(M,1);
% rng(1)
% % % %% normal distribution of etaSquared
% % % distChoice = 'normal';
% % % m_etaSquared = 6.45
% % % etaSquared = m_etaSquared + sigma_etaSquare*randn(M,1);
% % %% uniform distribution of etaSquared
% % distChoice = 'uniform';
% % etaSquared = a + (b-a)*rand(M,1);
% %% truncated normal distribution of etaSquared
% m_etaSquared = 6.45
% sigma = 1;
% pd = makedist('Normal','mu',m_etaSquared,'sigma',sigma);
% trunc_normal = truncate(pd,a,b);
% distChoice = 'truncNormal';
% etaSquared = random(trunc_normal,M,1);
% 

% %% for nomenclature of figures
% aval = strcat('a',num2str(a));
% aval = strrep(aval,'.','K');
% bval = strcat('b',num2str(b));
% bval = strrep(bval,'.','K');
% % only for truncated normal distribution needed
% sigmaval = strcat('sigma',num2str(sigma));
% sigmaval = strrep(sigmaval,'.','K');

%% for nomenclature of figures
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

%% Define gray scale for plots
mygray = flipud(gray((length(Q)+5)/2)); % if length(Q) uneven otherwise adapt grayscale
mycolors = [mygray(2:end-1,:) ;flipud(mygray(2:end-2,:))];


%% Initialize vector of fold points that end the attracting part of the critical manifold
firstCritEtaP2 = find(etaSquared_sort>3,1);
% %% comment if if firstCritEtaP2>N(1)
% numToDraw = N(length(Q))- N(1);
% uncomment if firstCritEtaP2>N(1)
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
% uncomment if firstCritEtaP2>N(1)
    if etaSquared_sort(j)>3
        %% Critical manifold splits into three parts: calculate fold point where attracting part ends
        x_crit(j-max(firstCritEtaP2,N(1))+1) = 2/3 +sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
        y_crit(j-max(firstCritEtaP2,N(1))+1) = h0(x_crit(j-max(firstCritEtaP2,N(1))+1),etaSquared_sort);
        
       %% strong fold point
       x_crit2(j-max(firstCritEtaP2,N(1))+1) = 2/3 -sqrt(4/9 - (1+etaSquared_sort(j))/(3*etaSquared_sort(j)));
        y_crit2(j-max(firstCritEtaP2,N(1))+1) = h0(x_crit2(j-max(firstCritEtaP2,N(1))+1),etaSquared_sort);   
    end
    % plot 3D probabilistic bifurcation curves
    plot3(y,etaSquared_sort(j)*ones(size(y)),x,'color',mycolors(i,:));
    view(315,45) % rotate for a better angle of view
end
end

h0 = @(x,etaSquared) x.*(1+3*(1-x).^2);
y=h0(x,etaSquared_sort);
plot3(y,3*ones(size(y)),x,'r','LineWidth',2);
c=colorbar;
mapColors = mycolors(3:end,:);
colormap(mapColors);
set(c, 'YTick', linspace(1/6, 5/6, 3));
c.TickLabels = num2cell(massLevels);

xlabel('Nondimensional freshwater flux');
ylabel('\eta^2');
zlabel('Salinity Difference');
ax.FontSize = 12;
ax.Interpreter = 'latex';
set(gca,'FontSize',12);

% % Save figure with 3D probabilistic representation of eta^2 dependent bifurcation curves in gray
% % scale
% savefig(fig1,strcat('probBifurcationCurves3D_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.fig'));
% saveas(fig1,strcat('probBifurcationCurves3D_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES','.eps'),'epsc');