function [nll,negScore,informationMatrix] = ...
         p_negLogLik(choices,iS,supportS,capPi,theta,delta,sp,sq,df,flowpayoffs,p_Bellman,p_fixedpoint,tolFixedPoint)

nSuppS = size(supportS,1);

[u0,u1] = flowpayoffs(supportS,theta,delta);
%scapPi = [capPi(:,1)+ sp*capPi(:,2) (1-sp-sq)*capPi(:,2) capPi(:,3)+sq*capPi(:,2)];
scapPi = [(1-sq)*capPi(:,1)+sp*capPi(:,2) sq*capPi(:,1)+(1-sp-sq)*capPi(:,2)+sp*capPi(:,3) sq*capPi(:,2)+(1-sp-sq)*capPi(:,3)+sp*capPi(:,4) sq*capPi(:,3)+(1-sp-sq)*capPi(:,4)+sp*capPi(:,5) sq*capPi(:,4)+(1-sp)*capPi(:,5)];

[capU0,capU1] = p_fixedpoint(u0,u1,scapPi,df,tolFixedPoint,p_Bellman,[],[]);
deltaU = capU1-capU0;
pExit = 1./(1+exp(deltaU));
laggedChoices = [zeros(1,size(choices,2));choices(1:end-1,:)];
p = choices + (1-2*choices).*pExit(iS+nSuppS*laggedChoices);
nll = -sum(sum(log(p)));

if nargout>=2

    d00 = df*scapPi*diag(pExit(:,1));
    d01	= df*scapPi*diag(pExit(:,2));
    d10	= df*scapPi-d00;
    d11	= df*scapPi-d01;
    dPsi_dUbar = [[d00;d00;zeros(2*nSuppS,nSuppS)] [zeros(2*nSuppS,nSuppS);d01;d01] ...
                  [d10;d10;zeros(2*nSuppS,nSuppS)] [zeros(2*nSuppS,nSuppS);d11;d11]];
	dPsi_dTheta = [[zeros(2*nSuppS,1);ones(2*nSuppS,1)] [zeros(2*nSuppS,1);supportS;supportS] [zeros(2*nSuppS,1);-ones(nSuppS,1);zeros(nSuppS,1)]];

    dUbar_dTheta   = (eye(4*nSuppS)-dPsi_dUbar)\dPsi_dTheta;
	dDeltaU_dTheta    = dUbar_dTheta(2*nSuppS+1:4*nSuppS,:)-dUbar_dTheta(1:2*nSuppS,:);	
    
    nTheta	= size(dUbar_dTheta,2);    
    negFirmScores = repmat((1-2*choices).*(1-p),[1 1 nTheta]);
    for i=1:nTheta
        negFirmScores(:,:,i) = negFirmScores(:,:,i).*dDeltaU_dTheta(iS+nSuppS*laggedChoices+2*(i-1)*nSuppS);
    end
    negFirmScores=squeeze(sum(negFirmScores,1));
    negScore = sum(negFirmScores)';
end

if nargout==3
    informationMatrix = zeros(nTheta,nTheta);
    for n=1:size(negFirmScores,1)
        informationMatrix = informationMatrix + negFirmScores(n,:)'*negFirmScores(n,:);
    end
end