% Donn?es bathym?triques provenant de
% http://www.gebco.net/data_and_products/gridded_bathymetry_data/
% taille de la grille : 30s d'arc, soit un demi-mille marin (926 m) en
% longitude

%% Section 1
% Chargement et affichage des donn?es altim?triques
clear all;
Z=load('map.asc');
[I J]=size(Z);                      % Map dimensions

%% section 2
% Initialisation des param?tres
N = 300;                 % Number of time steps.
t = 1:1:N;              % Time.
v0=1;                   % initial speed along x1
% x = zeros(N,2);         % Hidden states.
% y = zeros(N,1);         % Observations.
x(1,1) = 110;            % Initial state.
x(1,2) = 350;           % Initial state.
Rreal = 10^2;           % Measurement noise real variance.
R=5*Rreal;                % Measurement noise used for estimation
Qreal = [0.1 0;0 10];   % Process noise real variance.
Q = 2*Qreal;              % Process noise in y variance used for estimation
initVar = [100 0;0 100];    % Initial variance of the states.
numSamples=200;             % Number of Particles per time step.

%% Section 3
% Generation de la trajectoire et des mesures
for t=2:N
    x(t,:)=x(t-1,:)+[v0 0]+randn(1,2)*sqrt(Qreal); % trajectory (process) (time, position-> [x,y])
end
%filtrage de la trajectoire, pour "adoucissement"
alpha=0.01;
b=1-alpha;
a=[1 -alpha];
x=filter(b,a,x);    % removes noise from trajectory in y
%mesures
v = sqrt(Rreal)*randn(N,1); % measurement noise
for t=1:N
    y(t,1) = interp2(Z,x(t,1),x(t,2)) + v(t,1); % measurement (depth -> z)
end


%% Section 4 - Single execution of particle filter

display_it = true;         % Parameter to display/hide iterations in map
rsmpl_method='uniform';     % Resampling methdod choice

% Particules initiales (prior)
xxu=zeros(N,2,numSamples);
xu=sqrt(initVar)*randn(2,numSamples);
q=ones(1,numSamples);
xu(1,:)=xu(1,:)+x(1,1);     % Creation of 100 realizations of a gaussian random var following N(x(t=1,1), initVar) --> X dimension
xu(2,:)=xu(2,:)+x(1,2);     % Creation of 100 realizations of a gaussian random var following N(x(t=1,2), initVar) --> Y dimension


% Update et prediction
% clf(1);
% clf(2);
if ishandle(2), clf(2); end
if ishandle(3), clf(3); end
if ishandle(4), clf(4); end
if ishandle(3), clf(5); end
    
% --- Loop execution ---
it_loop;
% ----------------------

post_treatment;

%% Section 5 - Average performance tests
% This section performs the partoclar filter run several times, and gets
% the statistics out of them, for finding a mean performance in these runs
% -------------------------------------------------------------------------

% --- Settings ---
display_it = true;         % Parameter to display/hide iterations in map
rsmpl_method='uniform';     % Resampling methdod choice
M = 1;                      % Number of algorithm iterations

% Initialization of memory vectors
norms = zeros(M,4);         % Norm vector
n_eff_stat = zeros(M,1);    % Neff vector

for i=1:M
    disp("Ongoing iteration: "+i)
    
    % --- Initialization of particles (prior) ---
    xxu=zeros(N,2,numSamples);
    xu=sqrt(initVar)*randn(2,numSamples);
    q=ones(1,numSamples);
    xu(1,:)=xu(1,:)+x(1,1);     % Creation of 100 realizations of a gaussian random var following N(x(t=1,1), initVar) --> X dimension
    xu(2,:)=xu(2,:)+x(1,2);     % Creation of 100 realizations of a gaussian random var following N(x(t=1,2), initVar) --> Y dimension

    % --- Loop execution ---
    it_loop;
    
    % --- Image generation ---
    for fig=1:5
    if ishandle(fig), clf(fig); end     % Clears all figures
    end
    post_treatment;
    
    % --- Data saving ---
    x_data = x(:,1);
    y_data = x(:,2);
    x_data = x_data(1:end-1)';      % Reshaping trajectory X data
    y_data = y_data(1:end-1)';      % Reshaping trajectory Y data
    
    norms(i,1) = norm(EstX_hist(1,:)-x_data);  % X-axis difference 2-norm
    norms(i,2) = norm(EstX_hist(1,:)-x_data,'Inf');  % X-axis difference INF-norm
    norms(i,3) = norm(EstX_hist(2,:)-y_data);  % Y-axis difference 2-norm
    norms(i,4) = norm(EstX_hist(2,:)-y_data,'Inf');  % Y-axis difference INF-norm
    n_eff_stat(i,1) = norm(Neff_hist,'Inf');   % Effective particle inf norm
end

norms
mean_vec = mean(norms,1);
disp("Mean of norm-2 X-axis differences = " + mean_vec(1,1))
disp("Mean of norm-INF X-axis differences = " + mean_vec(1,2))
disp("Mean of norm-2 Y-axis differences = " + mean_vec(1,3))
disp("Mean of norm-INF Y-axis differences = " + mean_vec(1,4))
disp("Mean of maximal effective particle number = " + mean(n_eff_stat(isfinite(n_eff_stat))))