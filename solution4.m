clc;close all;clear
%% Condition 1
clc;close all;clear
a1=4.5; a2=1; a3=0.5; a4=0.4; a5=0.7; a6=0.2; a7=0.4;

% [(a1+sqrt(a1^2+4*a1))/a1;0;0]
x0 = [1;1;1];
t = 0:0.01:10000;
sol = simulate(t, x0, a1, a2, a3, a4, a5, a6, a7);
sol(end,1)
sol(end,2)
sol(end,3)

figure
for iPlot= 1:3
  subplot(1, 3, iPlot);
  plot(t, sol(:,iPlot),'LineWidth',2);
  legend("states vs time");
  title(["State", num2str(iPlot)])
end

% figure
% plot(sol(:,1), sol(:,2),'LineWidth',2);
% legend("x1 vs x2")
% 
% figure
% plot(sol(:,1), sol(:,3),'LineWidth',2);
% legend("x1 vs x3")
% 
% figure
% plot(sol(:,2), sol(:,3),'LineWidth',2);
% legend("x2 vs x3")

% figure
% plot3(sol(:,1), sol(:,2), sol(:,3),'LineWidth',2)
% xlabel("x1");ylabel("x2");zlabel("x3")
%% Condition 2
clc;close all;clear
a1=3.5; a2=1; a3=5; a4=0.4; a5=0.7; a6=0.1; a7=0.1;


x0 = [1;1;0];
t = 0:0.01:1000;
sol = simulate(t, x0, a1, a2, a3, a4, a5, a6, a7);

sol(end,1)
sol(end,2)
sol(end,3)

figure
for iPlot= 1:3
  subplot(1, 3, iPlot);
  plot(t, sol(:,iPlot),"LineWidth",2);
  title(["State", num2str(iPlot)])
end

% figure
% plot(sol(:,1), sol(:,2));
% legend("x1 vs x2")
% 
% figure
% plot(sol(:,1), sol(:,3));
% legend("x1 vs x3")
% 
% figure
% plot(sol(:,2), sol(:,3));
% legend("x2 vs x3")
% 
% figure
% plot3(sol(:,1), sol(:,2), sol(:,3))
% xlabel("x1")
% ylabel("x2")
% zlabel("x3")

%% Condition 3
clc;close all;clear
a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;

x0 = [1;0;1];
t = 0:0.01:1000;
sol = simulate(t, x0, a1, a2, a3, a4, a5, a6, a7);

sol(end,1)
sol(end,2)
sol(end,3)

figure
for iPlot= 1:3
  subplot(1, 3, iPlot);
  plot(t, sol(:,iPlot), 'LineWidth',2);
end

% figure
% plot(sol(:,1), sol(:,2));
% legend("x1 vs x2")
% 
% figure
% plot(sol(:,1), sol(:,3));
% legend("x1 vs x3")
% 
% figure
% plot(sol(:,2), sol(:,3));
% legend("x2 vs x3")
% 
% figure
% plot3(sol(:,1), sol(:,2), sol(:,3))
% xlabel("x1")
% ylabel("x2")
% zlabel("x3")


%% Function Definition

function sol=simulate(t, x0, a1, a2, a3, a4, a5, a6, a7)
    [t, sol] = ode45(@tumor,t, x0, [], a1, a2, a3, a4, a5, a6, a7);
end

function dydx = tumor(t, x, a1, a2, a3, a4, a5, a6, a7)
    dydx = [1+a1*x(1)*(1-x(1))-a2*x(1)*x(2);
            a3*x(2)*x(3)-a4*x(2);
            a5*x(3)*(1-x(3))-a6*x(2)*x(3)-a7*x(3)];
end