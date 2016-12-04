function piHat = estimatePi(iS,nSuppS)
nPeriods=size(iS,1);
%サンプルから遷移確率推定
for i=1:nSuppS
    for j=1:nSuppS
        piHat(i,j)=sum(sum((iS(2:nPeriods,:)==j)&(iS(1:nPeriods-1,:)==i)))/sum(sum((iS(1:nPeriods-1,:)==i)));
    end
end