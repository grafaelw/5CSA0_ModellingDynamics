% After playing around with the parameters, it can be concluded that there
% are 2 different fixed-point that are stable for H*. However, both of the
% M and R tend to converge to zero in this sense. We tried with using [1 1
% 1] as the initial value, and then raising it up until [10 10 10]. For the
% delay, we use from tau=10 up to tau=1000 and from that notion, we know
% that there are two stable fixed-points in the time-delayed dynamics.
%
%
clc,close all,clear all
tau=[1,10,100]; % Time delay regarding the equations

    
for iPlot=1:length(tau)
   sol = dde23(@ddetumor,[tau(iPlot), tau(iPlot), tau(iPlot)],ddetumorhist(j),[0 100]);
   figure(iPlot)
   plot(sol.x,sol.y,"LineWidth",2)
   title(strcat('Tau = ', num2str(tau(iPlot))));
%         xlabel('time t', 'FontWeight','bold', 'FontSize',12);
   ylabel('$\dot{x}(t)$', 'FontWeight','bold', 'FontSize',12, "Interpreter","latex");
   legend(["Tumor cells", "Hunting cells" ,"Resting cells"], "FontWeight","bold")
   set(gca, "FontSize", 14, "TitleFontWeight", "bold")
end



% ylim([0 15])


% --------------------------------------------------------------------------

function s = ddetumorhist(i)
% Testing using fixed point for constant history function.
%     s = [0.8263 1.7325 0.0833]';
      s = [1 1 1]';
end

% --------------------------------------------------------------------------

function dydt = ddetumor(~,y,Z)
% Differential equations function for DDEX1.
a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;
% Defining the lag function
ylag1 = Z(:,1);
ylag2 = Z(:,2);
ylag3 = Z(:,3);
dydt = [ 1+a1*y(1)*(1-y(1))-a2*y(1)*y(2)
        a3*ylag2(2)*ylag3(3)-a4*y(2)
        a5*y(3)*(1-y(3))-a6*y(2)*y(3)-a7*y(3)];
end

