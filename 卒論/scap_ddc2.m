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
%nSuppS = 3;
nSuppS = 5;
supportS = (1:nSuppS)';
%capPi = [0.9 0.1 0 ; 0.05 0.9 0.05 ; 0 0.1 0.9];
p = 0.1; 
q = 0.1;
%scapPi = [capPi(:,1)+ p*capPi(:,2) (1-p-q)*capPi(:,2) capPi(:,3)+q*capPi(:,2)]
capPi = [0.9 0.1 0 0 0 ; 0.05 0.9 0.05 0 0 ; 0 0.05 0.9 0.05 0 ; 0 0 0.05 0.9 0.05 ; 0 0 0 0.1 0.9];
scapPi = [(1-q)*capPi(:,1)+p*capPi(:,2) q*capPi(:,1)+(1-p-q)*capPi(:,2)+p*capPi(:,3) q*capPi(:,2)+(1-p-q)*capPi(:,3)+p*capPi(:,4) q*capPi(:,3)+(1-p-q)*capPi(:,4)+p*capPi(:,5) q*capPi(:,4)+(1-p)*capPi(:,5)]

theta = [-1;10];
delta = 1.0;
df = 0.95;	
[u0,u1] = flowpayoffs(supportS,theta,delta); 
[capU0,capU1] = p_fixedpoint(u0,u1,scapPi,df,tolFixedPoint,@p_Bellman,[],[]);
[ocapU0,ocapU1] = fixedPoint(u0,u1,capPi,df,tolFixedPoint,@Bellman,[],[]);
deltaU = capU1-capU0;
odeltaU = ocapU1-ocapU0;
%サンプル生成（モンテカルロシミュレーション）
for se=1:100
s = se
[choices,iS] = p_simulatedata(deltaU,capPi,nPeriods,nFirms,se);
%推定する関数の生成
objectiveFunction = @(parameters)negLogLik(choices,iS,supportS,scapPi,parameters(1:2),parameters(3),...
                                         df,@flowpayoffs,@Bellman,@fixedPoint,tolFixedPoint);
  %objectiveFunction = @(parameters)p_negLogLik(choices,iS,supportS,capPi,parameters(1:2),parameters(3),...
      %parameters(4),parameters(5),df,@flowpayoffs,@p_Bellman,@p_fixedpoint,tolFixedPoint);
for seed=1:9
rng(seed,'twister');  
sn = seed
%初期値
startvalues(:,1) = [-1.0;10.0;1.0];
startvalues(:,seed+1) = [-2 + 2*rand(1,1);7 + 6*rand(1,1);0 + 2*rand(1,1)];
%startvalues(:,1) = [-1.0;5.0;1.0;0.3;0.3];
%startvalues(:,seed+1) = [-2 + 2*rand(1,1);4 + 2*rand(1,1);2*rand(1,1);0.1 + 0.3*rand(1,1);0.1 + 0.2*rand(1,1)];
options = optimset('Display','off','TolFun',1E-3,'TolX',1E-3);
%fminuncでやったら推定値がおかしくなった
[maxLikEstimates,fval,exitflag] = fminunc(objectiveFunction,startvalues(:,seed),options);
%upp = [0 0 0 1 1;0 0 0 -1 -1;0 0 0 1 0;0 0 0 0 1;0 0 0 -1 0;0 0 0 0 -1];
%uppn = [1;0;1;1;0;0];
%[maxLikEstimates,fval,exitflag] = fmincon(objectiveFunction,startvalues(:,seed),upp,uppn,[],[],[],[],[],options);
[~,~,informationMatrix] = objectiveFunction(maxLikEstimates);
[nll,~,~]=objectiveFunction(maxLikEstimates);
pnll=0.00002*nll;
standardErrors =[ diag(sqrt(inv(informationMatrix)))];
MOMO(:,seed)=[maxLikEstimates;pnll];
[~,MOP]=min(MOMO(4,:));
MON=MOMO(:,MOP);
end

MO(:,se)=MON;
end
reportEstimates = mean(MO');
reportstd = std(MO');
%結果
disp('Summary of Results');
disp('--------------------------------------------');
disp('      true      estim      ste.');
disp([[theta;delta;0] reportEstimates' reportstd' ]);

%確率遷移の推定
%piHat = estimatePi(iS,nSuppS);