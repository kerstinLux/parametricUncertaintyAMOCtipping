function salDiff = uq_redStommelCessiODEModel(X,time)
% implementation of the reduced Stommel Cessi ODE model
% The implemented structure is taken from https://www.uqlab.com/inversion-predator-prey, last checked: 13.07.2021.

firstTime = time(1);
lastTime = time(end);    % Duration time of simulation.

nSteps = numel(time); % Number of timesteps
%% Initialize
nReal = size(X,1); 

% Parameters
etaSquared = X(:,1);
p = X(:,3);

% Initial conditions
initialSalDiff = X(:,2);

%% Solve equation

% solver options (smaller tolerance)
odeOpts = odeset('RelTol',1e-4,'AbsTol',1e-7);
%odeOpts = odeset('RelTol',2e-3,'AbsTol',2e-6);

% for loop to solve equations with multiple initial values and parameters
salDiff = zeros(nReal,nSteps);
for i = 1:nReal
    % setup diff equations 
    diffEq=@(t,x) p(i) - x(1)*(1+etaSquared(i)*(1-x(1))^2);
    % solve using numerical ODE solver 45
    [t,sol] = ode45(diffEq,[firstTime lastTime],initialSalDiff(i),odeOpts);
    % interpolate solution to specified timesteps
    interpSolSalDiff = interp1(t,sol(:,1),time);
    % assign solution
    salDiff(i,:) = [interpSolSalDiff'];
end

