function [ S_STRUCT , pChange, ChangeRate , Deflection ] = JMiPOD_fast (C_STRUCT , Payload)
% -------------------------------------------------------------------------
% J_MiPOD Embedding       |      February 2021       |      version 1.0
% -------------------------------------------------------------------------
% INPUT:
%  - C_STRUCT    - Struct representing JPEG compressed image (or path to JPEG file)
%  - Payload     - Embedding payload in bits per non-zeros AC DCT coefs (bpnzAC).
% OUTPUT:
%  - S_STRUCT    - Resulting stego jpeg STRUCT with embedded payload
%  - pChange     - Embedding change probabilities. 
%  - ChangeRate  - Average number of changed pixels
%  - Deflection  - Overall deflection. 
% -------------------------------------------------------------------------
% Copyright (c) 2020
% Remi Cogranne, UTT (Troyes University of Technology)
% All Rights Reserved.
% -------------------------------------------------------------------------
% This code is provided by the author under Creative Common License 
% (CC BY-NC-SA 4.0) which, as explained on this webpage
% https://creativecommons.org/licenses/by-nc-sa/4.0/
% Allows modification, redistribution, provided that:
% * You share your code under the same license ;
% * That you give credits to the authors ;
% * The code is used only for non-commercial purposes (which includes
% education and research)
% -------------------------------------------------------------------------
% 
% The authors hereby grant the use of the present code without fee, and 
% without a written agreement under compliance with aforementioned  
% and provided and the present copyright notice appears in all copies.
% The program is supplied "as is," without any accompanying services from
% the UTT or the authors. The UTT does not warrant the operation of the 
% program will be uninterrupted or error-free. The end-user understands 
% that the program was developed for research purposes and is advised not 
% to rely exclusively on the program for any reason. In no event shall 
% the UTT or the authors be liable to any party for any consequential 
% damages. The authors also fordid any practical use of this code for 
% communication by hiding data into JPEG images. 
% -------------------------------------------------------------------------
% Contact: remi.cogranne@utt.fr
%          Septembre 2020
% -------------------------------------------------------------------------
% References to be mentioned to give appropriate credits to the authors:
% R.Cogranne, Q.Giboulot & P.Bas "Efficient Steganography in JPEG Images by
% Minimizing Performance of Optimal Detector" under review
% ---------- %
% R.Cogranne, Q.Giboulot & P.Bas "Steganography by Minimizing Statistical
% Detectability: The cases of JPEG and Color Images", ACM IH&MMSec, 2020,
% Denver, CO, USA, pp. 161--167 https://doi.org/10.1145/3369412.3395075.
% ---------- %
% This work has been developped for :
% R.Cogranne, Q.Giboulot & P.Bas "ALASKA-2 : Challenging Academic Research
% on Steganalysis with Realistic Images", IEEE WIFS, 2020, NYC, USA 
% ---------- %
% Based on the prior work for spatial domain images :
% V. Sedighi, R. Cogranne, and J. Fridrich. "Content-Adaptive Steganography by Minimizing Statistical Detectability."
% IEEE Transactions on Information Forensics and Security, vol. 11, no. 2, 221--234, 2016 
% -------------------------------------------------------------------------
% Let us wrap it up by acknowledging the great inspiration we owe to Prof.
% Jessica Fridrich and to Vahid Sedighi with whom the spatial version of
% this method had been developed.
% Part of this code have been based on Matlab Implementation of MIPOD

% Read the JPEG image if needed
if ischar( C_STRUCT )
    C_STRUCT = jpeg_read( C_STRUCT );
end

%First let us decompress the image back into the spatial domain
[ C_SPATIAL , MatDCT ] = RCdecompressJPEG(C_STRUCT);

% Compute Variance in spatial domain ....
WienerResidual = C_SPATIAL - wiener2(C_SPATIAL,[2,2]);
Variance = VarianceEstimationDCT2D(WienerResidual,3,3);

% ... and apply the covariance transformation to DCT domain
% funVar = @(x) reshape( diag(MatDCT*diag(x(:))*MatDCT')  , 8 , 8 ) ./ ( C_STRUCT.quant_tables{1}.^2 );
% VarianceDCT = blkproc(Variance,[8 8],funVar);

% In this code we replaced the blkproc with nested loops and simplied covariance linear transformation
MatDCTq = MatDCT.^2;
Qvec = C_STRUCT.quant_tables{1}(:);
for idx=1:64 , MatDCTq(idx,:) = MatDCTq(idx,:)./ Qvec(idx).^2; end

VarianceDCT = zeros(size(C_SPATIAL));
for idxR=1:8:size( Variance, 1)
    for idxC=1:8:size( Variance, 2)
        x = Variance(idxR:idxR+7 , idxC:idxC+7);
        VarianceDCT(idxR:idxR+7 , idxC:idxC+7) = reshape( MatDCTq * x(:) , 8,8);
    end
end
VarianceDCT(VarianceDCT<1e-5) = 1e-5;

% Compute Fisher information and smooth it
FisherInformation = 1./VarianceDCT.^2;

%Post Filter
tmp = zeros(size( FisherInformation ) + 16);
tmp(9:end-8, 9:end-8) = FisherInformation;
tmp(1:8, :) = tmp(9:16, :);
tmp(end-7:end, :) = tmp(end-15:end-8, :);
tmp( : , 1:8, :) = tmp( : , 9:16);
tmp( : , end-7:end, :) = tmp( : , end-15:end-8);
FisherInformation =  tmp(1:end-16 , 1:end-16) + tmp(9:end-8 , 1:end-16) * 3 + tmp(17:end , 1:end-16) + tmp(1:end-16 , 9:end-8) * 3 + tmp(9:end-8 , 9:end-8) * 4 + tmp(17:end , 9:end-8) * 3 + tmp(1:end-16 , 17:end) + tmp(9:end-8 , 17:end) * 3 + tmp(17:end , 17:end) ; 

% Compute embedding change probabilities and execute embedding
FI = FisherInformation(:)';

C_COEFFS = C_STRUCT.coef_arrays{1};
S_COEFFS = C_COEFFS;

% Compute message nbnzAC and message lenght (in NATS !!)
nzAC = sum(C_COEFFS(:)~=0) - sum(sum(sum(C_COEFFS(1:8:end,1:8:end,:)~=0 ) ) );
messageLenght = round(Payload * nzAC * log(2));

[ beta ] = TernaryProbs(FI,messageLenght);

% Simulate embedding
beta = 2 * beta;
r = rand(1,numel(C_COEFFS));
ModifPM1 = (r < beta);                % Cover elements to be modified by +-1
r = rand(1,numel(VarianceDCT));
S_COEFFS(ModifPM1) = C_COEFFS(ModifPM1) + 2*(round(r(ModifPM1))) - 1; % Modifying X by +-1
S_COEFFS(S_COEFFS>1024) = 1024;                    % Taking care of boundary cases
S_COEFFS(S_COEFFS<-1023)   = -1023;
ChangeRate = sum(ModifPM1(:))/numel(C_COEFFS); % Computing the change rate
pChange = reshape(beta,size(C_COEFFS));

S_STRUCT = C_STRUCT;
S_STRUCT.coef_arrays{1} = S_COEFFS;
S_STRUCT.dc_huff_tables = {};
S_STRUCT.ac_huff_tables = {};  

Deflection = sum( pChange(:) .* FI(:) );

end

% Beginning of the supporting functions

%JPEG Image decompression and obtention of a single (64x64) matrix for DCT transform
function [ imDecompress , dct64_mtx ] = RCdecompressJPEG(imJPEG)
    %DCT matrix
    [cc,rr] = meshgrid(0:7);
    T = sqrt(2 / 8) * cos(pi * (2*cc + 1) .* rr / (2 * 8));
    T(1,:) = T(1,:) / sqrt(2);

    %DCT transformation as a single 64x64 matrix
    dct64_mtx = zeros(64,64);
    for i=1:64 ; dcttmp=zeros(8); dcttmp(i)=1; TTMP =  T*dcttmp*T'; dct64_mtx(:,i) = TTMP(:); end

    %Apply image decompression
    imDecompress = zeros( [ size( imJPEG.coef_arrays{1} ), numel( imJPEG.coef_arrays ) ] );
    DCTcoefs = imJPEG.coef_arrays{1};
    QM = imJPEG.quant_tables{1};
% Replace the blkproc use with nested loops
%     funIDCT = @(x) T'*(x.*QM)*T ;
%     imDecompress = blkproc(DCTcoefs,[8 8],funIDCT);
    for idxR=1:8:size( DCTcoefs, 1)
        for idxC=1:8:size( DCTcoefs, 2)
            imDecompress(idxR:idxR+7 , idxC:idxC+7) = T'*(DCTcoefs(idxR:idxR+7 , idxC:idxC+7).*QM)*T;
        end
    end
end

% Estimation of the pixels' variance based on a 2D-DCT (trigonometric polynomial) model
% Same function as in MiPOD
function EstimatedVariance = VarianceEstimationDCT2D(Image, BlockSize, Degree)
    % verifying the integrity of input arguments
    if ~mod(BlockSize,2)
        error('The block dimensions should be odd!!');
    end
    if (Degree > BlockSize)
        error('Number of basis vectors exceeds block dimension!!');
    end

    % number of parameters per block
    q = Degree*(Degree+1)/2;

    % Build G matirx
    BaseMat = zeros(BlockSize);BaseMat(1,1) = 1;
    G = zeros(BlockSize^2,q);
    k = 1;
    for xShift = 1 : Degree
        for yShift = 1 : (Degree - xShift + 1)
            G(:,k) = reshape(idct2(circshift(BaseMat,[xShift-1 yShift-1])),BlockSize^2,1);
            k=k+1;
        end
    end

    % Estimate the variance
    PadSize = floor(BlockSize/2*[1 1]);
    I2C = im2col(padarray(Image,PadSize,'symmetric'),BlockSize*[1 1]);
    PGorth = eye(BlockSize^2) - (G*((G'*G)\G'));
    EstimatedVariance = reshape(sum(( PGorth * I2C ).^2)/(BlockSize^2 - q),size(Image));
end

% Computing the embedding change probabilities
% Updated from MiPOD for speed purposes
function [beta] = TernaryProbs(FI, payload)
    load('ixlnx3_5.mat');

    % Initial search interval for lambda
    [L, R] = deal (0.1 , 10);

    fL = h_tern(1./invxlnx3_fast(L*FI,ixlnx3)) - payload;
    fR = h_tern(1./invxlnx3_fast(R*FI,ixlnx3)) - payload;
    % If the range [L,R] does not cover alpha enlarge the search interval
    while fL*fR > 0
        if fL > 0
            L = R;
            R = 2*R;
            fR = h_tern(1./invxlnx3_fast(R*FI,ixlnx3)) - payload;
        else
            R = L;
            L = L/2;
            fL = h_tern(1./invxlnx3_fast(L*FI,ixlnx3)) - payload;
        end
    end

    % Search for the labmda in the specified interval
    i = 0;
    M = (L+R)/2;
    fM = h_tern(1./invxlnx3_fast(M*FI,ixlnx3)) - payload;
    while (abs(fM)>max(2,payload/1000.0) && i<20)
        if fL*fM < 0,
            R = M; fR = fM;
        else
            L = M; fL = fM;
        end
        i = i + 1;
        M = (L+R)/2;
        fM = h_tern(1./invxlnx3_fast(M*FI,ixlnx3)) - payload;
    end
    % Compute beta using the found lambda
    beta = 1./invxlnx3_fast(M*FI,ixlnx3);
end


% Fast solver of y = x*log(x-2) paralellized over all pixels
% Updated from MiPOD for speed purposes
function x = invxlnx3_fast(y,f)
    i_large = y>1000;
    i_small = y<=1000;

    iyL = floor(y(i_small)/0.05)+1;
    iyR = iyL + 1;
    iyR(iyR>20000) = 20000;

    x = zeros(size(y));
    x(i_small) = f(iyL) + (y(i_small)-(iyL-1)*0.05).*(f(iyR)-f(iyL));

    z = y(i_large)./log(y(i_large)-2);
    z = y(i_large)./log(z-2);
    x(i_large) = y(i_large)./log(z-2);
end

% Ternary entropy function expressed in nats
function Ht = h_tern(Probs)
    p0 = 1-2*Probs;
    P = [p0(:);Probs(:);Probs(:)];
    H = -(P .* log(P));
    H((P<eps)) = 0;
    Ht = nansum(H);
end
