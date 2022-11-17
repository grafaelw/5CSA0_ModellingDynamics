%% Lyapunov function for Condition 3
clc;close all;clear

% Parameters
a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;

J = @(M, H, R)[-H*a2-a1*(2*M - 1),-M*a2,0;......
                0, R*a3-a4,H*a3;.......
                0,-R*a6, -a7-H*a6-R*a5-a5*(R - 1)];

% Fixed Points
Mstar = 0.8263; Hstar = 1.7325; Rstar = 0.0833;
P = lyap(J(Mstar, Hstar, Rstar),eye(3))
V = @(M, H, R) [M-Mstar, H-Hstar, R-Rstar]*P*[M-Mstar, H-Hstar, R-Rstar]';

Q=J(Mstar, Hstar, Rstar)'*P+P*J(Mstar, Hstar, Rstar)
eig(Q)


%% Function Definition

% function sol=simulate(t, x0, a1, a2, a3, a4, a5, a6, a7)
%     [t, sol] = ode45(@tumor,t, x0, [], a1, a2, a3, a4, a5, a6, a7);
% end
% 
% function dydx = tumor(t, x, a1, a2, a3, a4, a5, a6, a7)
%     dydx = [1+a1*x(1)*(1-x(1))-a2*x(1)*x(2);
%             a3*x(2)*x(3)-a4*x(2);
%             a5*x(3)*(1-x(3))-a6*x(2)*x(3)-a7*x(3)];
% end