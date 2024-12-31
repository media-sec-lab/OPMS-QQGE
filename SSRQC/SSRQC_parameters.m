function par = SSRQC_parameters()
% the iteration is fixed to 2 in this paper
% other parameters are using the dafaulting settings 

%% Parameters for optimization
par.IterNum = 2;      %Iteration is fixed to 2 in this paper
par.rho = 0.05;
par.delta = 0.5;      % 0.5 instead of 0.2

%% Parameter for group construction
par.BlockSize = 8;              % size of the block for DCT
par.PatchSize = 6;              % size of the patch for modelling
par.Profile = 'fast';

if strcmp(par.Profile, 'normal')
    par.SlidingDis   =  2;
elseif strcmp(par.Profile, 'fast')
    par.SlidingDis   =  5;
end
par.SearchWin = 20;
par.ArrayNo = 30;

% In the implementation, we do not explicitly set a value for lambda, but
% instead set the value for sqrt(2*lambda*K*rho)/N which is Factor as
% follows
par.Factor = par.PatchSize + sqrt(par.ArrayNo);

%% Other parameters
par.Qfactor = 0.4;


