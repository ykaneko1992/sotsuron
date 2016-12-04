function [choices,iS] = simulatedata(deltaU,capPi,nPeriods,nFirms)
%5*1s—ñ
nSuppX = size(capPi,1);
oneMinPi = eye(nSuppX)-capPi';
pInf     = [oneMinPi(1:nSuppX-1,:);ones(1,nSuppX)]\[zeros(nSuppX-1,1);1];
% 
iS = randomDiscrete(pInf*ones(1,nFirms));
%‹É’l•ª•z‚É]‚¤ƒÃ¶¬
deltaEpsilon = random('ev',zeros(1,nFirms),ones(1,nFirms))-random('ev',zeros(1,nFirms),ones(1,nFirms));
%ŠÖŒW‰‰ŽZŽq‚Í^‹U‚Å0,1‚ð•Ô‚·
choices  = deltaU(iS,1)' > deltaEpsilon;

for t = 2:nPeriods
	iS = [iS;randomDiscrete(capPi(iS(end,:),:)')]; 
	deltaEpsilon  = random('ev',zeros(1,nFirms),ones(1,nFirms))-random('ev',zeros(1,nFirms),ones(1,nFirms));
	choices = [choices;(deltaU(iS(end,:)+nSuppX*choices(end,:)) > deltaEpsilon)];
end