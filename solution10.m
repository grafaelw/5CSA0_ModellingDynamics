%% Linearisation of the system around desired and fixed points
clc;close all;clear

clc;close all;clear
a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;

Mdesired = 0.01; Hdesired = 102.97; Rdesired = 1/12; udesired = 288527/18000; % solved by Matlab

A = @(M, H, R)[-H*a2-a1*(2*M - 1), -M*a2, 0;
                0, R*a3-a4, H*a3;
                0, -R*a6, -a7-H*a6-R*a5-a5*(R - 1)];
B = [0;0;1];

sysDesign = ss(A(Mdesired,Hdesired,Rdesired),B,[1,1,1],0);
sysFP = ss(A(0.8263,1.7325,0.0833),B,[1,1,1],0);
K = [1 10 1];
sys1 = ss(A(Mdesired,Hdesired,Rdesired)-(B*K),zeros(3,1),[1 1 1],0)

pole(sysDesign)

% Checking for controllability
Co = ctrb(sysDesign);
rank(Co);

% Need to find: A+BK

%% Simulation
clc;close all;clear
a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;
x0 = [100000;0.001;1;0];
t = 0:0.01:10000;
sol = simulate(t, x0, a1, a2, a3, a4, a5, a6, a7);

figure
plot(t,sol,'LineWidth',2)
%         xlabel('time t', 'FontWeight','bold', 'FontSize',12);
ylabel('Density', 'FontWeight','bold', 'FontSize',12);
legend(["Tumor cells", "Hunting cells" ,"Resting cells", "input"], "FontWeight","bold")
set(gca, "FontSize", 14, "TitleFontWeight", "bold")

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
    [t, sol] = ode45(@tumor,t,x0, [], a1, a2, a3, a4, a5, a6, a7);
end

function [dydx,u] = tumor(t, x, a1, a2, a3, a4, a5, a6, a7)
    K = [1 1 1];
    Mstar = 0.01; Hstar = 102.97; Rstar = 1/12; ustar = 288527/18000;
    u = ustar-K*([x(1)-Mstar;x(2)-Hstar;x(3)-Rstar]);
    dydx = [1+a1*x(1)*(1-x(1))-a2*x(1)*x(2);
            a3*x(2)*x(3)-a4*x(2);
            a5*x(3)*(1-x(3))-a6*x(2)*x(3)-a7*x(3)+u;
            u];
end