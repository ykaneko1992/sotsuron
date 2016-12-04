function y = randomDiscrete(p)
nSupp = size(p,1);
nVar = size(p,2);
uniformDraws = ones(nSupp-1,1)*random('unif',zeros(1,nVar),ones(1,nVar));
cumulativeP = cumsum(p);
y = sum([ones(1,nVar);cumulativeP(1:nSupp-1,:)<=uniformDraws]);
