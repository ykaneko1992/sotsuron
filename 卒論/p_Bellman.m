function [capU0,capU1] = p_Bellman(capU0,capU1,u0,u1,scapPi,df)
% U_a(S_t,a_{t-1})=u_a(S_t,a_{t-1})+df*E[R_a(s)|S_t]
% R_a(s)=log{{exp{U_0(s,a)}+exp{U_1(s,a)}}
r0 = log(exp(capU0(:,1))+exp(capU1(:,1)));
r1 = log(exp(capU0(:,2))+exp(capU1(:,2)));

capU0 = [u0 u0] + df*scapPi*r0*[1 1];
capU1 = u1 + df*scapPi*r1*[1 1];