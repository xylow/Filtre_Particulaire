%% Post-treatment script - Plot figures (and update data recording)

if ishandle(1)
    figure(1)
    contour(Z,50); 
    hx=line(x(:,1),x(:,2),'LineStyle','-','LineWidth',2,'Marker','none','Color','black');
    hestx=line(EstX_hist(1,:),EstX_hist(2,:),'LineStyle','-','LineWidth',2,'Marker','none','Color','red');
    xlim([min(EstX_hist(1,1),x(1,1)), max(EstX_hist(1,end),x(1,end))])
    ylim([min(min(min(EstX_hist(2,:),x(:,2)))) max(max(max(EstX_hist(2,:),x(:,2))))])
end

figure(3)
plot(1:N-1,Neff_hist)
title("Effective number of particles through time")
xlabel("time [s]")
ylabel("Neff")

figure(4)
plot(1:N-1,EstX_hist(1,:),'color','r'); hold on;
plot(1:N,x(:,1),'color','b')
legend('estimation','trajectory')
title("Estimation in the X axis through time")
xlabel("time [s]")
ylabel("X")

figure(5)
plot(1:N-1,EstX_hist(2,:),'color','r'); hold on;
plot(1:N,x(:,2),'color','b')
legend('estimation','trajectory')
title("Estimation in the Y axis through time")
xlabel("time [s]")
ylabel("Y")