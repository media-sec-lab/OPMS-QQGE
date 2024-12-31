function  [x] = SSRQC_Solver_Deblocking(Inoi, par)

y = Inoi;                   % noisy image y
x = Inoi;                   % initilize the denoised image x

% load parameters
sigma_e = par.nSig;         % sigma_e
rho = par.rho;              % rho
delta = par.delta;          % delta
IterNum = par.IterNum;

% augxilary variables
b = zeros(size(x));         % SBI variable b
D_a = x;                    % D_a: D * alpha, dictionary reconstructed image

for i  =  1 : IterNum 
   %% Step1: Solve Eq. (23)
    w = D_a - b;                % the equation above Eq. (27)
    
    % Calculate beta
    % Eq. (29), variance of the noise s
    vd = sigma_e^2-(mean(mean((w - y).^2)));      % mean: ||w-y||_2^2/N^2; typo in the paper
    if (i == 1)
        sigma_s2 = abs(vd);
    else
        sigma_s2 = abs(vd) * delta^2;
    end    
    % Eq. (28), updated beta
    beta = rho*sigma_s2; 
    
    % Eq. (26), update the denoised image x
    x = (beta * y + sigma_e^2*(D_a + b))/(beta + sigma_e^2);
    
    %% Step 2
    r = x - b;                  % the equation below Eq. (30)   
    % For recalculating beta
    vd = sigma_e^2-(mean(mean((r - y).^2)));    
    if (i == 1)
        sigma_s  = sqrt(abs(vd));
    else
        sigma_s  = sqrt(abs(vd)) * delta;
    end     
       
    % Eq. (31), min 1/2beta ||D_a - r||_2^2 + lambda ||alpha||_0
    D_a0 = SSRQC_Dict_learning(r, par, sigma_s);
    
    % Eq. (32), min 1/2 ||D_a - D_a0||_2^2 + phi(D_a)
    D_a = SSRQC_BDCT_project_onto_QCS(D_a0, par);
    
    %% Step 3: Eq. (25)
    b = b + (D_a - x);
    
end
