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
%capPi = [0.8 0.2 0 0 0 ; 0.2 0.6 0.2 0 0 ; 0 0.2 0.6 0.2 0 ; 0 0 0.2 0.6 0.2 ; 0 0 0 0.2 0.8];
%scapPi = [0.85 0.1 0.05 0 0 ; 0.1 0.4 0.4 0.1 0 ; 0.1 0.2 0.4 0.2 0.1 ; 0 0.1 0.4 0.4 0.1; 0 0 0.05 0.1 0.85];
theta = [-1;1];
delta = 3.0;
df = 0.95;	
[u0,u1] = flowpayoffs(supportS,theta,delta); 
[capU0,capU1] = p_fixedpoint(u0,u1,scapPi,df,tolFixedPoint,@p_Bellman,[],[]);
[ocapU0,ocapU1] = fixedPoint(u0,u1,capPi,df,tolFixedPoint,@Bellman,[],[]);
deltaU = capU1-capU0;
odeltaU = ocapU1-ocapU0;
%サンプル生成（モンテカルロシミュレーション）
[choices,iS] = simulatedata(deltaU,capPi,nPeriods,nFirms);
%推定する関数の生成
objectiveFunction = @(parameters)negLogLik(choices,iS,supportS,scapPi,parameters(1:2),parameters(3),...
                                         df,@flowpayoffs,@Bellman,@fixedPoint,tolFixedPoint);
  %objectiveFunction = @(parameters)p_negLogLik(choices,iS,supportS,capPi,parameters(1:2),parameters(3),...
      %parameters(4),parameters(5),df,@flowpayoffs,@p_Bellman,@p_fixedpoint,tolFixedPoint);  
%初期値
%startvalues = [-1.0;10.0;1.0];
%startvalues = [-2 + (2)*rand(1,1);5 + (15-5)*rand(1,1);2*rand(1,1)]
startvalues = [-1.5;15.0;2.0;0.1;0.4];
options = optimset('Display','iter','TolFun',1E-3,'TolX',1E-3);
%fminuncでやったら推定値がおかしくなった
[maxLikEstimates,fval,exitflag] = fminunc(objectiveFunction,startvalues,options);
[~,~,informationMatrix] = objectiveFunction(maxLikEstimates);
[nll,~,~]=objectiveFunction(maxLikEstimates);
pnll=0.00002*nll
standardErrors =[ diag(sqrt(inv(informationMatrix)));0;0];
%結果
disp('Summary of Results');
disp('--------------------------------------------');
disp('      true     start     estim      ste.');
disp([[theta;delta;p;q] startvalues maxLikEstimates standardErrors]);

%確率遷移の推定
%piHat = estimatePi(iS,nSuppS);