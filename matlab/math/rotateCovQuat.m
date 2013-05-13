function [m,cov] = rotateCovQuat(mu, Sigma, q)

R = (quat2rm(q));

m = (R*mu')';
cov = R*Sigma*R';

%DEBUG:
% x.pos = [0 0 0];
% x.rot = q;
% m_q = getXfrmPointQuat(x, mu);
% 
% m
% m_q