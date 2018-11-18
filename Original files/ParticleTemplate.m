% Données bathymétriques provenant de
% http://www.gebco.net/data_and_products/gridded_bathymetry_data/
% taille de la grille : 30s d'arc, soit un demi-mille marin (926 m) en
% longitude

%% Section 1
% Chargement et affichage des données altimétriques
clear all;
Z=load('map.asc');
[I J]=size(Z);
figure(1);contour(Z,50);
set(gca,'Nextplot','replace');

%% section 2
% Initialisation des paramètres
N = 500;                % Number of time steps.
t = 1:1:N;              % Time.
v0=1;                   % speed along x1
% x = zeros(N,2);         % Hidden states.
% y = zeros(N,1);         % Observations.
x(1,1) = 30;            % Initial state.
x(1,2) = 170;           % Initial state.
Rreal = 10^2;           % Measurement noise real variance.
R=Rreal;             % Measurement noise used for estimation
Qreal = [0.1 0;0 10];   % Process noise real variance.
Q = Qreal;                 % Process noise variance used for estimation
initVar = [0 0;0 0]; % Initial variance of the states.
numSamples=100;         % Number of Particles per time step.

%% Section 3
% Génération de la trajectoire et des mesures
for t=2:N,
    x(t,:)=x(t-1,:)+[v0 0]+randn(1,2)*sqrt(Qreal); % trajectory (process)
end;
%filtrage de la trajectoire, pour "adoucissement"
alpha=0.01;
b=1-alpha;
a=[1 -alpha];
x=filter(b,a,x);
%mesures
v = sqrt(Rreal)*randn(N,1); % measurement noise
for t=1:N,
    y(t,1) = interp2(Z,x(t,1),x(t,2)) + v(t,1); % measurement
end;


%% section 4
% Particules initiales (prior)
xxu=zeros(N,2,numSamples);
xu=sqrt(initVar)*randn(2,numSamples);
q=ones(1,numSamples);
xu(1,:)=xu(1,:)+x(1,1);
xu(2,:)=xu(2,:)+x(1,2);

hx=line(x(:,1),x(:,2),'LineStyle','-','Marker','none','Color','black');%tracé de la trajectoire complète
hxpos=line(0,0,'LineStyle','none',...
    'Marker','o','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','black');%tracé de la positioncourante
hxu=line(0,0,'LineStyle','none','Marker','o');%tracé des particules
hxell=line(0,0,'LineStyle','-','Marker','none');%tracé de l'ellipse


%% Section 5
% Update et prédiction
for t=1:N-1
    %Predict
    %from the set of particles xu generate a new set xu
    xu=xu;
    
    %Importance wheights
    for k=1:numSamples
        m(k)=interp2(Z,xu(1,k),xu(2,k));%mesures prédites pour chaque particules
    end;
    %from the set of weights q compute the new set of weights q
    q=q; %REMPLACER !
    [ii jj]=find(xu(1,:)>J | xu(1,:)<1 | xu(2,:)>I | xu(2,:)<1 );
    q(jj)=0; %Elimine les éventuelles particules "hors du cadre"
    q=q./sum(q);
    
    %Resampling
    Neff=1/sum(q.^2)
    if Neff<0.75*numSamples
        %Resamplpling
        method='none';
        switch method
            case 'none'
                xur=xu;
                q=q;
            case 'multinomial'
                xur=xu;%REMPLACER !!
                q=ones(1,numSamples)/numSamples;
        xu=xur;
    end
    
    %Stockage dans xxu
    xxu(t,:,:)=xu;
    %mise a jour affichages
    pause(0.001);
    set(hxu,'Xdata',reshape(xu(1,:),1,numSamples),...
        'Ydata',reshape(xu(2,:),1,numSamples));
    %     set(hx,'Xdata',x(1:t,1),'Ydata',x(1:t,2));
    set(hxpos,'Xdata',x(t,1),'Ydata',x(t,2));
    %Ellipsoid
    X0=mean(xu,2);
    M=(xu-X0*ones(1,length(xu)))*(xu-X0*ones(1,length(xu)))'/length(xu);
    [U,S,V]=svd(M);
    a=0:0.1:2*pi;
    ell0=[sqrt(S(1,1))*cos(a);sqrt(S(2,2))*sin(a)];
    ell=X0*ones(1,length(a))+V'*ell0;
    set(hxell,'Xdata',ell(1,:),'Ydata',ell(2,:),'LineWidth',2,'Color','r'); 
end

