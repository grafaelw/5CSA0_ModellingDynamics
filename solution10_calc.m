clear; clc
syms y z u

a1=3; a2=1; a3=4.8; a4=0.4; a5=3.7; a6=1.9; a7=0.1;
x=0.01;

eqn1 = 0==1+a1*x*(1-x)-a2*x*y;
eqn2 = 0==a3*y*z-a4*y;
eqn3 = 0==a5*z*(1-z)-a6*y*z-a7*z+u;

eqns = [eqn1, eqn2, eqn3];
S = solve(eqns,[y z, u],'ReturnConditions',true);
S.y
S.z
S.u