%% Condition 1 ==> Asymptotically stable
clc;close all;clear
a1=4.5; a2=1; a3=0.5; a4=0.4; a5=0.7; a6=0.2; a7=0.4;
J = @(M, H, R)[-H*a2-a1*(2*M - 1),-M*a2,0;
                0, R*a3-a4,H*a3;
     0,-R*a6, -a7-H*a6-R*a5-a5*(R - 1)];

% Fixed Points
Mstar = 1.1870; Hstar = 0; Rstar = 0.4286;
J(Mstar, Hstar, Rstar)
eig(J(Mstar, Hstar, Rstar))

% Since we know that the stability of a fixed point occurs at the local minima, 
% we can try to do a linearisation around the fixed points by taking the Jacobian 
% matrix as it is the first derivative of a function. Negative eigen-values
% correspond to it implies that the matrix is negative definite and it is
% indeed the 

%% Condition 2 ==> Unstable
clc;close all;clear
a1=3.5; a2=1; a3=5; a4=0.4; a5=0.7; a6=0.1; a7=0.1;
J = @(M, H, R)[-H*a2-a1*(2*M - 1),-M*a2,0;
                0, R*a3-a4,H*a3;
     0,-R*a6, -a7-H*a6-R*a5-a5*(R - 1)];

% Fixed Points
Mstar1 = 1.2318; Hstar1 = 0; Rstar1 = 0.8571;
J(Mstar1, Hstar1, Rstar1)
eig(J(Mstar1, Hstar1, Rstar1))


%% Condition 3 ==> Asymptotically stable
clc;close all;clear
a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;
J = @(M, H, R)[-H*a2-a1*(2*M - 1),-M*a2,0;
                0, R*a3-a4,H*a3;
     0,-R*a6, -a7-H*a6-R*a5-a5*(R - 1)];

% Fixed Points
Mstar2 = 0.8263; Hstar2 = 1.7325; Rstar2 = 0.0833;
J(Mstar2, Hstar2, Rstar2)
eig(J(Mstar2, Hstar2, Rstar2))







