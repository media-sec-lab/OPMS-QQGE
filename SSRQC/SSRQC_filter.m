function IF=SSRQC_filter(I,qtab)
% SSRQC filter with fixed iteration being 2

par = SSRQC_parameters();          % Set parameters
par.nim = I;                       % load the JPEG noisy image y
par.QTable = qtab;        % quantization table for 8*8 block
par.C_q  = blkproc(par.nim , [8, 8], 'dct2');    % apply DCT to each block of y      
meanQuant=mean(mean(par.QTable(1:3,1:3)));
par.nSig = sqrt(0.69*meanQuant^1.3);             % Gaussian variance for the quantization noise model 
IF = SSRQC_Solver_Deblocking(par.nim, par);

end