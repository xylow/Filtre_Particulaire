Neff_hist = zeros(1,N-1);
EstX_hist = zeros(2,N-1);

% Initialize graphical objects
if display_it
    % Level curves plotting
    figure(1);contour(Z,50);            
    set(gca,'Nextplot','replace');
    % Complete trajectory tracing
    hx=line(x(:,1),x(:,2),'LineStyle','-','LineWidth',2,'Marker','none','Color','black');
    % Current position tracing
    hxpos=line(0,0,'LineStyle','none',...
        'Marker','o','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','black');
    % Particle tracing
    hxu=line(0,0,'LineStyle','none','Marker','o');     
    % Uncertainty ellipse tracing
    hxell=line(0,0,'LineStyle','-','Marker','none');   
end

for t=1:N-1
    % --- Predict ---
    % Propagate the set of particles xu according to covariance matrix Q
    xu=xu + [v0 0]' + sqrt(Q)*randn(2,numSamples);      % Random dispersion 
    
    %Importance wheights
    for k=1:numSamples                      % k-th particle's
        m(k)=interp2(Z,xu(1,k),xu(2,k));    % Interpolated predicted measure
    end
    
    % From the set of weights q propagate the new set of weights q
    % For this we use a Gaussian distribution: N ~ (measure y , R matrix)
    q=q.*exp(-1/(2*R)*(m-y(t,1).*ones(size(m))).^2)/sqrt(2*pi*R);
    % Finds particles out of the map zone
    [ii,jj]=find(xu(1,:)>J | xu(1,:)<1 | xu(2,:)>I | xu(2,:)<1 );
    q(jj)=0;        % Elimine les eventuelles particules "hors du cadre"
    q=q./sum(q);    % Weight normalization
    
    %Resampling
    Neff=1/sum(q.^2);       % Calculates effective particle number
    Neff_hist(t) = Neff;    % Saves to memory
    
    % --- Resampling ---
    if Neff<0.75*numSamples % If we have a low number of particles
        switch rsmpl_method
            case 'none'
                xur=xu;
                q=q;
            case 'uniform'
                [xur,q] = resample(xu,q); % Gets new particles and weights
            case 'multinomial'
                
        end
                xu=xur;
    end
    
    %Stockage dans xxu
    xxu(t,:,:)=xu;
    %mise a jour affichages
    
    %Ellipsoid
    X0=mean(xu,2);
    EstX_hist(:,t) = X0;
    M=(xu-X0*ones(1,length(xu)))*(xu-X0*ones(1,length(xu)))'/length(xu);
    [U,S,V]=svd(M);
    a=0:0.1:2*pi;
    ell0=[sqrt(S(1,1))*cos(a);sqrt(S(2,2))*sin(a)];
    ell=X0*ones(1,length(a))+V'*ell0;
    if display_it
        pause(0.001);
        set(hxu,'Xdata',reshape(xu(1,:),1,numSamples),...
            'Ydata',reshape(xu(2,:),1,numSamples));
        %     set(hx,'Xdata',x(1:t,1),'Ydata',x(1:t,2));
        set(hxpos,'Xdata',x(t,1),'Ydata',x(t,2));
        % Ellipse
        set(hxell,'Xdata',ell(1,:),'Ydata',ell(2,:),'LineWidth',2,'Color','r');
%         figure(2);
%         plot(1:numSamples, sort(q));
%         histogram(sort(q))
%         title("Particle weight in iteration "+t)
    end
end