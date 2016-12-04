function [capU0,capU1] = p_fixedpoint(u0,u1,scapPi,df,tolFixedPoint,p_Bellman,capU0,capU1)

% initial parameter 
nSuppS= size(scapPi,1);
if isempty(capU0) 
	capU0 = zeros(nSuppS,2);
end
if isempty(capU1)
	capU1 = zeros(nSuppS,2);
end

%solving the fixed point of the Bellman eq
inU0 = capU0+2*tolFixedPoint;
inU1 = capU1+2*tolFixedPoint;
% コメント貰う　最高iteration回数設定してforで回したほうがいい？
while (max(max(abs(inU0-capU0)))>tolFixedPoint) || (max(max(abs(inU1-capU1)))>tolFixedPoint);
	inU0 = capU0;
	inU1 = capU1;
	[capU0,capU1] = p_Bellman(inU0,inU1,u0,u1,scapPi,df);
end