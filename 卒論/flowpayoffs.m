%define utility functions
function [u0,u1] = flowpayoffs(supportS,theta,delta)
nSuppS = size(supportS,1);
% u0 = 0 if a_t = 0
% u1 = theta_FC - theta_ec(1-a_i{t-1})+ theta_RS*S_t if a_t = 1
u0 = zeros(nSuppS,1);
% 15”{‚ªŒÀ“x
u1 = ([ones(nSuppS,1) supportS]*[theta(1); delta]*[1 1]-theta(2)*ones(nSuppS,1)*[1 0]);
