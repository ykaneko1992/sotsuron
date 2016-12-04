clear
for seed=1:2
rng(seed,'twister');  
sn = seed
%‰Šú’l
%startvalues(:,1) = [-1.0;5.0;1.0];
%startvalues(:,seed+1) = [-2 + 2*rand(1,1);5 + 10*rand(1,1);0 + 2*rand(1,1)];
startvalues(:,1) = [-1.0;5.0;1.0;0.3;0.3];
startvalues(:,seed+1) = [-2 + 2*rand(1,1);4 + 2*rand(1,1);2*rand(1,1);0.2 + 0.2*rand(1,1);0.2 + 0.2*rand(1,1)]
end