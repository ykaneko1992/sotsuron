clear
%{
nperiods
nFirms
tolFixedpoint
nSuppS=number of market state
supportS=value of market state
capPi=transition probability
theta =[(fixed effect);(Entry cost)]
delta = theta_RS
df = discount factor
	%}
nPeriods = 500;
nFirms = 100;
tolFixedPoint = 1e-100;
nSuppS = 3;
%nSuppS = 5;
supportS = (1:nSuppS)';
capPi = [0.9 0.1 0 ; 0.05 0.9 0.05 ; 0 0.1 0.9];
p = 0.3; 
q = 0.3;
scapPi = [capPi(:,1)+ p*capPi(:,2) (1-p-q)*capPi(:,2) capPi(:,3)+q*capPi(:,2)];
%capPi = [0.9 0.1 0 0 0 ; 0.05 0.9 0.05 0 0 ; 0 0.05 0.9 0.05 0 ; 0 0 0.05 0.9 0.05 ; 0 0 0 0.1 0.9];
%scapPi = [(1-q)*capPi(:,1)+p*capPi(:,2) q*capPi(:,1)+(1-p-q)*capPi(:,2)+p*capPi(:,3) q*capPi(:,2)+(1-p-q)*capPi(:,3)+p*capPi(:,4) q*capPi(:,3)+(1-p-q)*capPi(:,4)+p*capPi(:,5) q*capPi(:,4)+(1-p)*capPi(:,5)];

%theta = [-1;1];
%delta = 1.0;
theta=[-0.85;5];
delta = 0.82;
df = 0.95;	
[u0,u1] = flowpayoffs(supportS,theta,delta); 
[capU0,capU1] = p_fixedpoint(u0,u1,scapPi,df,tolFixedPoint,@p_Bellman,[],[])
[ocapU0,ocapU1] = fixedPoint(u0,u1,capPi,df,tolFixedPoint,@Bellman,[],[]);
deltaU = capU1-capU0;
odeltaU = ocapU1-ocapU0;
%サンプル生成（モンテカルロシミュレーション）
for se=1:100
s = se;
[choices,iS] = p_simulatedata(deltaU,capPi,nPeriods,nFirms,se);
size(choices);
k=find(choices);
a=size(k);
AB(se) = a(1);
end
shit=mean(AB)/50000