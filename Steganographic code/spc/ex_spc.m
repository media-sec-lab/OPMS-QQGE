% clc; clear;

% parameter setting
% n = 2^15 + 2^14;
% n = 2^20 + 1;
n = 2^20;      % cover length
alpha = 0.5;   % relative payload
L = 1;         % list size - default is 1 - drives the complexity/quality tradeof
wetcost = 1e10;   % wet cost for forbidding illegal modification

% generate random cover and message
% rand('seed',1234);
cover = int32(rand(1, n)*255);
m = round(n * alpha);           % number of message bits
message = uint8(rand(1, m));

% generate distortion profile (square)
rho = 1:n;
rho = rho .* rho;  
rho = rho(randperm(n)) / sum(rho);  % create a random permutation and normalize 
profile = single(rho);

% embed message
ary = 0;
while (ary ~= 2 && ary ~= 3 && ary ~= 4 && ary ~= 5)
    prompt = 'Please input \"2\" for binary or \"3\" for ternary or \"4\" for quaternary or \"5\" for pentary embedding.\n';
    ary = input(prompt);
end
% ary=2;
if ary==2      % binary embedding
    t1 = clock;
%     [dist, stego, n_msg_bits] = spc_lsb_pls_embed(cover, profile, message, L); % embed message
    [dist, stego, n_msg_bits] = spc_ml_pls_embed(cover, profile, message, L); % another interface for embedding message    
    t2 = clock;
    extr_msg = spc_ml_extract(stego, n_msg_bits);
elseif ary==3       % ternary embedding
    costs = zeros(3, n, 'single'); 
    costs(1,:) = profile;      % cost of changing by -1
    costs(3,:) = profile;      % cost of changing by +1
    % prevent invalid stego pixels
    costs(1,cover==0) = wetcost;
    costs(3,cover==255) = wetcost;   
    t1 = clock;
%     [dist, stego, n_msg_bits] = spc_pm1_pls_embed(cover, costs, message, L); % embed message
    [dist, stego, n_msg_bits] = spc_ml_pls_embed(cover, costs, message, L); % another interface for embedding message  
    t2 = clock;
    extr_msg = spc_ml_extract(stego, n_msg_bits);   
elseif ary==4       % quaternary embedding
    costs = zeros(4, n, 'single'); 
    costs(1,:) = profile;      % cost of changing by -1
    costs(3,:) = profile;      % cost of changing by +1
    costs(4,:) = profile;      % cost of changing by +2/-2
    % prevent invalid stego pixels
    costs(1,cover==0) = wetcost;
    costs(3,cover==255) = wetcost;   
    costs(4,cover==255) = wetcost;
    t1 = clock;
%     [dist, stego, n_msg_bits] = spc_pm1_pls_embed(cover, costs, message, L); % embed message
    [dist, stego, n_msg_bits] = spc_ary4_pls_embed(cover, costs, message, L); % another interface for embedding message  
    t2 = clock;
    extr_msg = spc_ml_extract(stego, n_msg_bits);
elseif ary==5        % pentary embedding
    costs = zeros(5, n, 'single'); 
    costs(1,:) = 4 * profile;      % cost of changing by -2
    costs(2,:) = profile;          % cost of changing by -1    
    costs(4,:) = profile;          % cost of changing by +1
    costs(5,:) = 4 * profile;      % cost of changing by +2    
    % prevent invalid stego pixels
    costs(1,cover==1) = wetcost;
    costs(2,cover==0) = wetcost;    
    costs(4,cover==254) = wetcost; 
    costs(5,cover==255) = wetcost;    
    t1 = clock;
%     [dist, stego, n_msg_bits] = spc_pm2_pls_embed(cover, costs, message, L); % embed message
    [dist, stego, n_msg_bits] = spc_ml_pls_embed(cover, costs, message, L); % another interface for embedding message   
    t2 = clock;
    extr_msg = spc_ml_extract(stego, n_msg_bits); 
end
interval = etime(t2,t1);
fprintf('\nTotal distortion = %.4f amd its embedding efficiency = %.4f.\nEmbedding time = %.4f sec and its throughput = %.4f Kbits/sec.\n', dist, m/dist, interval, n/interval/1024);
if all(extr_msg==message)
    fprintf('Message is embedded and extracted correctly.\n');
else
    error('Some errors occured in the extraction process.');
end
