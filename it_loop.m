Neff_hist = zeros(1,N-1);
EstX_hist = zeros(2,N-1);

for t=1:N-1
    %Predict
    %from the set of particles xu generate a new set xu
    xu=xu + [v0 0]' + sqrt(Q)*randn(2,numSamples);      % Random dispersion 
    
    %Importance wheights
    for k=1:numSamples  %k-th particle
        m(k)=interp2(Z,xu(1,k),xu(2,k));    %mesures predites pour chaque particules
    end
    %from the set of weights q compute the new set of weights q
    q=q.*exp(-1/(2*R)*(m-y(t,1).*ones(size(m))).^2)/sqrt(2*pi*R);
    [ii jj]=find(xu(1,:)>J | xu(1,:)<1 | xu(2,:)>I | xu(2,:)<1 );
    q(jj)=0; %Elimine les eventuelles particules "hors du cadre"
    q=q./sum(q);
    
    %Resampling
    Neff=1/sum(q.^2);
    Neff_hist(t) = Neff;
    if Neff<0.75*numSamples %|| isnan(Neff)
        %Resamplpling
        method='uniform';
        switch method
            case 'none'
                xur=xu;
                q=q;
            case 'uniform'
                [xur,q] = resample(xu,q);
            case 'multinomial'
                
        end
                xu=xur;
    end
    
    %Stockage dans xxu
    xxu(t,:,:)=xu;
    %mise a jour affichages
%     pause(0.001);
%     set(hxu,'Xdata',reshape(xu(1,:),1,numSamples),...
%         'Ydata',reshape(xu(2,:),1,numSamples));
%     %     set(hx,'Xdata',x(1:t,1),'Ydata',x(1:t,2));
%     set(hxpos,'Xdata',x(t,1),'Ydata',x(t,2));
    %Ellipsoid
    X0=mean(xu,2);
    EstX_hist(:,t) = X0;
    M=(xu-X0*ones(1,length(xu)))*(xu-X0*ones(1,length(xu)))'/length(xu);
    [U,S,V]=svd(M);
    a=0:0.1:2*pi;
    ell0=[sqrt(S(1,1))*cos(a);sqrt(S(2,2))*sin(a)];
    ell=X0*ones(1,length(a))+V'*ell0;
    %set(hxell,'Xdata',ell(1,:),'Ydata',ell(2,:),'LineWidth',2,'Color','r'); 
    
%     figure(2);
%     plot(1:numSamples, sort(q));
%     title("Particle weight in iteration "+t)
    
    sum(q);
end