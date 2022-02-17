% In this example, the reduced Stommel-Cessi model (Cessi, 1994) is
% calibrated against a time series of synthetically generated salinity data by a forward solve of the ODE and with a normal random perturbation added to it.
%
% The structure is taken from https://www.uqlab.com/inversion-predator-prey, last checked: 13.07.2021.
%
% The code is based on UQLab manual
    % P.-R. Wagner, J. Nagel, S. Marelli, B. Sudret, UQLab user manual – Bayesian inference for model calibration and
    % inverse problems, Report # UQLab-V1.3-113, Chair of Risk, Safety and Uncertainty Quantification, ETH Zurich,
    % Switzerland, 2019
% and predator prey example (https://www.uqlab.com/inversion-predator-prey, last checked: 13.07.2021).

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

%% 2 - FORWARD MODEL
%
% The forward model used for calibration is the solution of the 
% Stommel Cessi model reduced to the center manifold
% 
% The forward model computes the salinity difference for
% 51 large time instances, based on a forward ODE simulation with a small step
% size discretizing [0,T]
% The model takes as input parameters:
%
% $\eta^2$: ratio between advective and diffusive time scale
% $x_0$: fixed initial salinity difference
% $\mu = p$: nondimensional freshwater flux
%
% The computation is carried out by the function uq_redStommelCessiODEModel implemented based on the template |uq_predatorPreyModel|
% supplied with UQLab. For every set of input parameters, the function
% returns the salinity difference evolution in [0,T].

%%
% create column vector with discretization of synthetic data points from [0,T]
t0 = 0;
T=5;
time = (t0:0.1:T)';

%%
% Specify the forward models as a UQLab MODEL object:
ModelOpts.mHandle = @(x) uq_redStommelCessiODEModel(x,time); % ATTENTION: value of p fixed to PriorOpts.Marginals(3).Parameters as used for synthetic data
ModelOpts.isVectorized = true;

myForwardModel = uq_createModel(ModelOpts);

%% 3 - PRIOR DISTRIBUTION OF THE MODEL PARAMETERS
%
% To encode the available information about the model parameters before any experimental observations,
% a distribution is put on the parameter \eta^2 as follows:
%
% # $etaSquared \sim \mathcal{U}(0.6, 12.3)$
%
% Specify this prior distribution and the constant initial salinity difference
% as well as the nondimensional freshwater flux as a UQLab INPUT object:

PriorOpts.Marginals(1).Name = ('etaSquared');
%% Gaussian prior
% PriorOpts.Marginals(1).Type = 'Gaussian';
% PriorOpts.Marginals(1).Parameters = [4 1]; % second entry for Gaussian parameters is std
% %% Triangular prior
% PriorOpts.Marginals(1).Type = 'Triangular';
% PriorOpts.Marginals(1).Parameters = [0.6 12.3 6.45];
%% Uniform prior
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [0.6 12.3];
% PriorOpts.Marginals(1).Parameters = [0 5];

PriorOpts.Marginals(2).Name = ('initSalDiff');
PriorOpts.Marginals(2).Type = 'Constant';
PriorOpts.Marginals(2).Parameters = [0.4]; % close to attracting lower branch
% PriorOpts.Marginals(2).Parameters = [0.3]; % very close to equilibrium 0f 0.2729 at attracting lower branch
% PriorOpts.Marginals(2).Parameters = 1.05;
% PriorOpts.Marginals(2).Parameters = 0.5; % close to equilibrium value of ~0.3933
% PriorOpts.Marginals(2).Parameters = 1.4;

PriorOpts.Marginals(3).Name = ('p');
PriorOpts.Marginals(3).Type = 'Constant';
PriorOpts.Marginals(3).Parameters = 0.85;

myPriorDist = uq_createInput(PriorOpts);

%% 4 - MEASUREMENT DATA
%
%% synthetic data generation from ODE simulation
% etaSquared = 7.5; % Cessi (1994) value
etaSquared = 4;
p = PriorOpts.Marginals(3).Parameters;
% %% Calculate equilibrium salinity difference corresponding to chosen
% freshwater forcing value mu=p
% x0 = fzero(@(x) p-x.*(1+etaSquared*(1-x).^2),0);

% discretization
dt = 10^(-3);
tspan = t0:dt:T;
% initial value
x0 = PriorOpts.Marginals(2).Parameters;
% RHS
odefun = @(t,x) p - x(1)*(1+etaSquared*(1-x(1))^2);
% solve ODE
[~,x] = ode45(odefun,tspan,x0);
% noiseIntensity = 0.05;
noiseIntensity = 0.3;
% noiseIntensity = 0.5;
% noiseIntensity = 1;
xi = [0; noiseIntensity*randn(length(time)-1,1)]; % noise on data
synData = x(1:100:end)+ xi;
myData(1).y = synData';
myData(1).Name = 'Sal Diff data';
myData(1).MOMap = 1:T*10+1; % Output ID


%% 5 - DISCREPANCY MODEL
%
%% Give opportunity to learn variance sigma^2 of Gaussian discrepancy
% To infer the discrepancy variance, a lognormal prior is put on the
% discrepancy parameter:
%
% Specify this distribution in UQLab separately as an INPUT object:
SigmaOpts.Marginals(1).Name = 'SigmaP2';
SigmaOpts.Marginals(1).Type = 'Lognormal';
SigmaOpts.Marginals(1).Parameters = [-1 1];
% SigmaOpts.Marginals(1).Type = 'Uniform';
% SigmaOpts.Marginals(1).Parameters = [0 1];

SigmaDist = uq_createInput(SigmaOpts);

%%
% Assign this distribution to the discrepancy model options:
DiscrepancyOpts(1).Type = 'Gaussian';
DiscrepancyOpts(1).Prior = SigmaDist;

% %% Explicitly specify Gaussian discrepancy with fixed variance sigma^2
% DiscrepancyOpts(1).Type = 'Gaussian';
% DiscrepancyOpts(1).Parameters = 1e-4; % corresponds to sigma = 0.01

%% 6 - BAYESIAN ANALYSIS
%
%% 6.1 MCMC solver options
%
% To sample from the posterior distribution, the affine invariant ensemble
% algorithm is chosen, with $400$ iterations and $100$ parallel chains:
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.Steps = 400;
% Solver.MCMC.Steps = 1000;
Solver.MCMC.NChains = 100;
% Solver.MCMC.Sampler = 'MH';
% Solver.MCMC.Steps = 800;
% Solver.MCMC.NChains = 500;

%%
% Enable progress visualization during iteration for the etaSquared parameter
% Update the plots every $40$ iterations:
Solver.MCMC.Visualize.Parameters = [1];
Solver.MCMC.Visualize.Interval = 40;

%% 6.2 Posterior sample generation
%
% The options of the Bayesian analysis are gathered within a single
% structure with fields: 
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian model';
BayesOpts.Prior = myPriorDist;
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;

%%
% Perform and store in UQLab the Bayesian inversion analysis:
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

%%
% Print out a report of the results:
uq_print(myBayesianAnalysis)

%% 6.3 Posterior sample post-processing
%
% Diagnose the quality of the results,
% create a trace plot of the first parameter:
uq_display(myBayesianAnalysis, 'trace', 1)

%%
% Additionally, draw a sample of size $10^3$ from the prior and the prior and posterior
% predictive distributions: 
uq_postProcessInversion(myBayesianAnalysis,...
                        'prior', 1000,...
                        'priorPredictive', 1000,...
                        'posteriorPredictive', 1000);

%%
% *Note*: sampling prior predictive samples requires new
% model evaluations

%%
% Display the post processed results:
uq_display(myBayesianAnalysis)

%% Save sample from posterior distribution in .mat-file
PostSample3D = myBayesianAnalysis.Results.PostProc.PostSample;
% according to explanation in
% https://uqworld.org/t/looking-for-documentation-example-of-using-uqlab-results-of-bayesian-inference-for-forward-uq-or-further-bayesian-inference/690,
% last checked: 04.08.2021
PostSample2D = reshape(permute(myBayesianAnalysis.Results.PostProc.PostSample, [2 1 3]), size(PostSample3D, 2), []);
UQpostSample = PostSample2D(1,:);
UQpostStd = sqrt(myBayesianAnalysis.Results.PostProc.Percentiles.Var(1)); % posterior std
UQpostMean =  myBayesianAnalysis.Results.PostProc.Percentiles.Mean(1); % posterior mean
UQpostPercentiles = myBayesianAnalysis.Results.PostProc.Percentiles.Values(:,1);
% median
uq_postProcessInversion(myBayesianAnalysis,'percentiles', [0.5]);
UQpostMedian = myBayesianAnalysis.Results.PostProc.Percentiles.Values(1,1); % posterior median

%% Compute MAP of posterior distribution
uq_postProcessInversion(myBayesianAnalysis,'pointEstimate', 'MAP')
UQpostMap = myBayesianAnalysis.Results.PostProc.PointEstimate.X(1);

% logLikeliEval = myBayesianAnalysis.Results.LogLikeliEval;

% %% Save mat-file with posterior samples
% save('UQpostSample_y0_0K4_T5_etaP2_4_p0K85_noise0K3_priorUni_0K6_12K3_AIES_steps400_NChains100.mat','etaSquared','UQpostSample','UQpostStd','UQpostMean','UQpostMap','UQpostPercentiles');

