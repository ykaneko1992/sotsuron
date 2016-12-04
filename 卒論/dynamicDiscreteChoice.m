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
nPeriods = 100;
nFirms = 1000;
tolFixedPoint = 1e-10;
nSuppS = 5;
supportS = (1:nSuppS)';
capPi = [0.8 0.2 0 0 0 ; 0.2 0.6 0.2 0 0 ; 0 0.2 0.6 0.2 0 ;0 0 0.2 0.6 0.2; 0 0 0 0.2 0.8];
theta = [-1.5;-1.0]
delta = 1
df = 0.95	
[u0,u1] = flowpayoffs(supportS,theta,delta); 
[capU0,capU1] = fixedPoint(u0,u1,capPi,df,tolFixedPoint,@Bellman,[],[]);
deltaU = capU1-capU0;
%�T���v�������i�����e�J�����V�~�����[�V�����j
[choices,iS] = simulatedata(deltaU,capPi,nPeriods,nFirms);
%���肷��֐��̐���
objectiveFunction = @(parameters)negLogLik(choices,iS,supportS,capPi,parameters(1:2),parameters(3),...
                                           df,@flowpayoffs,@Bellman,@fixedPoint,tolFixedPoint)
%�����l
startvalues = [-1;0.5;0.5];
options = optimset('Display','iter','TolFun',1E-10,'TolX',1E-10);
%fminunc�ł�����琄��l�����������Ȃ���
[maxLikEstimates,fval,exitflag] = fminsearch(objectiveFunction,startvalues,options);
[~,~,informationMatrix] = objectiveFunction(maxLikEstimates);
standardErrors = diag(sqrt(inv(informationMatrix)));
%���ʕ\��
disp('Summary of Results');
disp('--------------------------------------------');
disp('      true     start     estim      ste.');
disp([[theta;delta] startvalues maxLikEstimates standardErrors]);
%�m���J�ڂ̐���
piHat = estimatePi(iS,nSuppS)