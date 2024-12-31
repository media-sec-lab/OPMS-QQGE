function  OPMS_QQGE_main(cover_img,stego_img,QF,payload,basecost,ary_num,coding_method)
%Side-information Estimated (SIE) JPEG Steganography
%OPMS-QQGE TIFS2023
%% Parameters:
% cover_img: path of cover image
% stego_img: path of stego image
% QF: quality factor
% payload: embedding payload (in bpnzAC)
% basecost: basic cost function (e.g. UNIWARD, JMiPOD)
% ary_num: k-ary embedding (e.g., 3 - ternary, 4 - quaternary, 5 - pentary)
% coding_method: coding method for stego generation (e.g., simulator, STC, SPC)

%Execution Example:
%OPMS_QQGE_main('cover-75\3.jpg','stego\3.jpg',75,0.5,'JMiPOD',4,'STC');

%% OPMS with 6 deblockers (AVG/Wie/SSRQC/DnCNN/DRUNet/SwinIR)
% generate modifications
img = jpeg_read(cover_img);
dct1 = img.coef_arrays{1};
[k,l]=size(dct1);
qtab = img.quant_tables{1}; 
q_matrix = repmat(qtab,[k/8 l/8]);
I_jpg=double(imread(cover_img));
rho=J_UNIWARD_COST_binary_mat(I_jpg,qtab);
wetConst = 10^13;     %wet cost
rhoP1_e=rho;
rhoM1_e=rho;
rhoP1_e(dct1>1023) = wetConst; 
rhoM1_e(dct1<-1023) = wetConst;
cover=double(dct1);
nzAC = nnz(cover)-nnz(cover(1:8:end,1:8:end));
N=numel(cover);
mlen=floor(nzAC*payload);
cover_1D=reshape(cover,N,1);
rhoP1_1D=reshape(rhoP1_e,N,1);
rhoM1_1D=reshape(rhoM1_e,N,1);
rho0=zeros(size(rhoP1_1D));
cost=[rhoP1_1D,rhoM1_1D,rho0];
cost=cost';
lambda=calc_lambda_curve_forward_lt1(cost,mlen,N);
pChangeP1 = (exp(-lambda .* rhoP1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
pChangeM1 = (exp(-lambda .* rhoM1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
randChange = rand(size(cover_1D));
stego_1D = cover_1D;
stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
stego=reshape(stego_1D,size(cover));
change=stego-dct1;

%% generate POSS
%Note: Please replace with the corresponding paths of the images deblocked by DRUNet and SwinIR
de_num=6;
de_path=strings(de_num,1);
dist_GFR=zeros(1,de_num);
if QF==75
    de_path(5)='Deblock2023/BOWS-256-QF75-DRUNet/';
    de_path(6)='Deblock2023/BOWS-256-QF75-SwinIR/';
elseif QF==95
    de_path(5)='Deblock2023/BOWS-256-QF95-DRUNet/';
    de_path(6)='Deblock2023/BOWS-256-QF95-SwinIR/';
elseif QF==50
    de_path(5)='Deblock2023/BOWS-256-QF50-DRUNet/';
    de_path(6)='Deblock2023/BOWS-256-QF50-SwinIR/';
end
%%
F=[1/9,1/9,1/9;1/9,1/9,1/9;1/9,1/9,1/9];     %3x3 average filter
IF=imfilter(I_jpg,F,'symmetric', 'conv', 'same');
IF_avg=IF;
IF_dct = bdct(round(IF)-128);
IF_dct_nr=IF_dct./q_matrix;
est_e=IF_dct_nr-dct1;
sgn_ee=sign(est_e);     %polarity of estimated rounding error
change2=change;
change2(change~=0 & sgn_ee~=0) = sgn_ee(change~=0 & sgn_ee~=0);
stego=dct1+change2;
S_STRUCT = img;
S_STRUCT.coef_arrays{1}=stego;
jpeg_write(S_STRUCT,stego_img);
dist_GFR(1) = calc_dist_feature(cover_img,stego_img,QF,'GFR');
%%
IF=wiener2(I_jpg,[3 3]);     %3x3 Wiener filter
IF_wie=IF;
IF_dct = bdct(round(IF)-128);
IF_dct_nr=IF_dct./q_matrix;
est_e=IF_dct_nr-dct1;
sgn_ee=sign(est_e);     %polarity of estimated rounding error
change2=change;
change2(change~=0 & sgn_ee~=0) = sgn_ee(change~=0 & sgn_ee~=0);
stego=dct1+change2;
S_STRUCT = img;
S_STRUCT.coef_arrays{1}=stego;
jpeg_write(S_STRUCT,stego_img);
dist_GFR(2) = calc_dist_feature(cover_img,stego_img,QF,'GFR');
%%
for i=3:6
    if i==3
%%     
        IF = SSRQC_filter(I_jpg,qtab);     %SSRQC filter
        IF_ssrqc=IF;
    elseif i==4
        if (QF==50 || QF==75 || QF==95)
            load(['DnCNN_Deblocker_QF',num2str(QF),'.mat']);
            IF = double(denoiseImage(uint8(I_jpg),net));
            IF_dncnn=IF;
        else
            error('Only applicable to QF50, QF75 and QF95.');
        end
    else
        path_test=de_path(i);
        imgs_dir = dir(cover_img);
        new_path=[char(path_test),imgs_dir.name(1:end-4),'.pgm'];
        IF= double(imread(new_path));
        if i==5
            IF_drunet=IF;
        elseif i==6
            IF_swinir=IF;
        end
    end
    IF_dct = bdct(round(IF)-128);
    IF_dct_nr=IF_dct./q_matrix;
    est_e=IF_dct_nr-dct1;
    sgn_ee=sign(est_e);     %polarity of estimated rounding error
    change2=change;
    change2(change~=0 & sgn_ee~=0) = sgn_ee(change~=0 & sgn_ee~=0);
    stego=dct1+change2;
    S_STRUCT = img;
    S_STRUCT.coef_arrays{1}=stego;
    jpeg_write(S_STRUCT,stego_img);
    dist_GFR(i) = calc_dist_feature(cover_img,stego_img,QF,'GFR');
%%
end
%% search for OPM
[val_gfr,col_gfr]=sort(dist_GFR);
col1_gfr=col_gfr(1);
if col1_gfr==1
    IF=IF_avg;
elseif col1_gfr==2
    IF=IF_wie;
elseif col1_gfr==3
    IF=IF_ssrqc;
elseif col1_gfr==4
    IF=IF_dncnn;
elseif col1_gfr==5
    IF=IF_drunet;
elseif col1_gfr==6
    IF=IF_swinir;
end
%%

%% embedding wity QQGE
img = jpeg_read(cover_img);
dct1 = img.coef_arrays{1};
[k,l]=size(dct1);
qtab = img.quant_tables{1}; 
q_matrix = repmat(qtab,[k/8 l/8]);
I_jpg=double(imread(cover_img));
IF_dct = bdct(round(IF)-128); 
IF_dct_nr=IF_dct./q_matrix;
est_e=IF_dct_nr-dct1;

if strcmp(coding_method,'STC')==1
    para=10;
elseif strcmp(coding_method,'SPC')==1
    para=1;
%     para=4;
end
% calculate basic costs
if strcmp(basecost,'UNIWARD')==1
    element=J_UNIWARD_COST_binary_mat(I_jpg,qtab);
    type=1;  %for cost
    if QF==75 || QF==50
        gamma=0.14+0.1*payload;
    elseif QF==95
        gamma=0.33-0.1*payload;
    end    
elseif strcmp(basecost,'JMiPOD')==1
    element = JMiPOD_fast_var(img);   
    type=2;  %for variance
    if QF==75 || QF==50
        gamma=0.09+0.1*payload;
    elseif QF==95
        gamma=0.28-0.1*payload;
    end
end
cover=double(dct1);
nzAC = nnz(cover)-nnz(cover(1:8:end,1:8:end));
N=numel(cover);
mlen=floor(nzAC*payload);
mlen_nat=round(mlen/log2(exp(1)));
est_e_p=gamma*sign(est_e);
[lambda,prob,mark]=calc_lambda_curve_forward_lt1_qgm_SI(element,mlen_nat,N,ary_num,type,est_e_p);
%% TQGE or QQGE or PQGE
if ary_num==3
    rhoP1 = log(prob(2,:)./prob(3,:));
    rhoM1 = log(prob(2,:)./prob(1,:));
    rhoP1_max=max(rhoP1(rhoP1~=inf));
    rhoM1_max=max(rhoM1(rhoM1~=inf));    
    rho_max=max(rhoP1_max,rhoM1_max);
    wetConst=max(rho_max*10^4,10^8);  %set wet cost
    rhoP1(isnan(rhoP1)) = wetConst;    
    rhoM1(isnan(rhoM1)) = wetConst;
    rhoP1(rhoP1==inf) = wetConst;    
    rhoM1(rhoM1==inf) = wetConst;   
    %embed
    cover_1D=reshape(cover,N,1);
    rhoP1_1D=reshape(rhoP1,N,1);
    rhoM1_1D=reshape(rhoM1,N,1);
    rho0=zeros(size(rhoP1_1D));
    cost=[rhoP1_1D,rhoM1_1D,rho0];
    cost=cost';    
    if strcmp(coding_method,'Sim')==1   %simulator        
        lambda=calc_lambda_curve_forward_lt1(cost,mlen,N);
        pChangeP1 = (exp(-lambda .* rhoP1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
        pChangeM1 = (exp(-lambda .* rhoM1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
        randChange = rand(size(cover_1D));
        stego_1D = cover_1D;
        stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
        stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
        stego=reshape(stego_1D,size(cover));
    elseif strcmp(coding_method,'STC')==1
        cover_stc=int32(cover_1D');
        cost_stc=single([rhoM1_1D';rho0';rhoP1_1D']);
        msg=randi([0,1],mlen,1);
        msg_stc=uint8(msg');
        try        
            [dist, stego_stc, n_msg_bits, l] = stc_pm1_pls_embed(cover_stc, cost_stc, msg_stc, para);
            stego_1D=double(stego_stc);
            stego=reshape(stego_1D,size(cover));
%             extr_msg = stc_ml_extract(int32(stego_stc), n_msg_bits, para);
%             success=0;
%             if all(extr_msg==msg_stc(1:numel(extr_msg)))
%                 success=1;
% %                 if numel(extr_msg)~=numel(msg_stc)
% %                     success=(numel(msg_stc)-numel(extr_msg))*100;
% %                 end
%             end
        catch
%             disp(cover_img);
%             delete(cover_img);
%             return;
            mark_success=0;
            lambda=calc_lambda_curve_forward_lt1(cost,mlen,N);
            %     lambda=calc_lambda(cost,mlen,N);
            pChangeP1 = (exp(-lambda .* rhoP1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
            pChangeM1 = (exp(-lambda .* rhoM1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
            randChange = rand(size(cover_1D));
            stego_1D = cover_1D;
            stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
            stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
            stego=reshape(stego_1D,size(cover));
        end
    elseif strcmp(coding_method,'SPC')==1
        cover_stc=int32(cover_1D');
        cost_stc=single([rhoM1_1D';rho0';rhoP1_1D']);
        msg=randi([0,1],mlen,1);
        msg_stc=uint8(msg');
        try
            [dist, stego_stc, n_msg_bits] = spc_pm1_pls_embed(cover_stc, cost_stc, msg_stc, para);
            stego_1D=double(stego_stc);
            stego=reshape(stego_1D,size(cover));
        catch
            mark_success=0;
            lambda=calc_lambda_curve_forward_lt1(cost,mlen,N);
            pChangeP1 = (exp(-lambda .* rhoP1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
            pChangeM1 = (exp(-lambda .* rhoM1_1D))./(1 + exp(-lambda .* rhoP1_1D) + exp(-lambda .* rhoM1_1D));
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
            randChange = rand(size(cover_1D));
            stego_1D = cover_1D;
            stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
            stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
            stego=reshape(stego_1D,size(cover));            
        end
    else 
        error('Without such coding method!');
    end   
elseif ary_num==4
    rhoP1 = log(prob(2,:)./prob(3,:));
    rhoM1 = log(prob(2,:)./prob(1,:));
    rhoPM2 = log(prob(2,:)./prob(4,:));    
    rho_max=max([rhoPM2(rhoPM2~=inf),rhoM1(rhoM1~=inf),rhoP1(rhoP1~=inf)]);
    wetConst=max(rho_max*10^4,10^8);  %set wet cost
    rhoP1(isnan(rhoP1)) = wetConst;    
    rhoM1(isnan(rhoM1)) = wetConst;
    rhoPM2(isnan(rhoPM2)) = wetConst; 
    rhoP1(rhoP1==inf) = wetConst;    
    rhoM1(rhoM1==inf) = wetConst;
    rhoPM2(rhoPM2==inf) = wetConst; 
    %embed
    cover_1D=reshape(cover,N,1);
    mark_1D=reshape(mark,N,1);
    rhoP1_1D=reshape(rhoP1,N,1);
    rhoM1_1D=reshape(rhoM1,N,1);
    rhoPM2_1D=reshape(rhoPM2,N,1);
    rho0=zeros(size(rhoP1_1D));
    cost=[rhoP1_1D,rhoM1_1D,rho0,rhoPM2_1D];
    cost=cost';    
    if strcmp(coding_method,'Sim')==1   %simulator
        [lambda,pp]=calc_lambda_curve_forward_lt1(cost,mlen,N);
        pChangeP1 = pp(1,:)';
        pChangeM1 = pp(2,:)';
        pChangePM2 = pp(4,:)';
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
        randChange = rand(size(cover_1D));
        stego_1D = cover_1D;
        stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
        stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
        stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==1) = stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==1) + 2;
        stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==0) = stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==0) - 2;
        stego=reshape(stego_1D,size(cover));
    elseif strcmp(coding_method,'STC')==1
        cover_stc=int32(cover_1D');
        rhoPM2_1D_2=rhoPM2_1D;
        rhoPM2_1D_2(mark==0)=-rhoPM2_1D_2(mark==0);   %cost>0 means +2 while cost<0 means -2 in stc_4ary code
        cost_stc=single([rhoM1_1D';rho0';rhoP1_1D';rhoPM2_1D_2']);
        msg=randi([0,1],mlen,1);
        msg_stc=uint8(msg');
        try
            [dist, stego_stc, n_msg_bits, l] = stc_ary4_pls_embed(cover_stc, cost_stc, msg_stc, para);
            stego_1D=double(stego_stc);
            stego=reshape(stego_1D,size(cover));
%             extr_msg = stc_ml_extract(int32(stego_stc), n_msg_bits, para);
%             success=0;
%             if all(extr_msg==msg_stc(1:numel(extr_msg)))
%                 success=1;
% %                 if numel(extr_msg)~=numel(msg_stc)
% %                     success=(numel(msg_stc)-numel(extr_msg))*100;
% %                 end
%             end            
        catch
%             disp(cover_img);
%             delete(cover_img);
%             return;
            mark_success=0;
            [lambda,pp]=calc_lambda_curve_forward_lt1(cost,mlen,N);
            pChangeP1 = pp(1,:)';
            pChangeM1 = pp(2,:)';
            pChangePM2 = pp(4,:)';
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
            randChange = rand(size(cover_1D));
            stego_1D = cover_1D;
            stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
            stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
            stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==1) = stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==1) + 2;
            stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==0) = stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==0) - 2;
            stego=reshape(stego_1D,size(cover));
        end
    elseif strcmp(coding_method,'SPC')==1
        cover_stc=int32(cover_1D');
        rhoPM2_1D_2=rhoPM2_1D;
        rhoPM2_1D_2(mark==0)=-rhoPM2_1D_2(mark==0);   %cost>0 means +2 while cost<0 means -2 in stc_4ary code
        cost_stc=single([rhoM1_1D';rho0';rhoP1_1D';rhoPM2_1D_2']);
        msg=randi([0,1],mlen,1);
        msg_stc=uint8(msg');
        try
            [dist, stego_stc, n_msg_bits] = spc_ary4_pls_embed(cover_stc, cost_stc, msg_stc, para);
            stego_1D=double(stego_stc);
            stego=reshape(stego_1D,size(cover));
        catch
            mark_success=0;
            [lambda,pp]=calc_lambda_curve_forward_lt1(cost,mlen,N);
            pChangeP1 = pp(1,:)';
            pChangeM1 = pp(2,:)';
            pChangePM2 = pp(4,:)';
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
            randChange = rand(size(cover_1D));
            stego_1D = cover_1D;
            stego_1D(randChange < pChangeP1) = stego_1D(randChange < pChangeP1) + 1;
            stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
            stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==1) = stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==1) + 2;
            stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==0) = stego_1D(randChange >= pChangeP1+pChangeM1 & randChange < pChangeP1+pChangeM1+pChangePM2 & mark_1D==0) - 2;
            stego=reshape(stego_1D,size(cover));            
        end
    else 
        error('Without such coding method!');
    end 
elseif ary_num==5
    rhoM2 = log(prob(3,:)./prob(1,:));
    rhoM1 = log(prob(3,:)./prob(2,:));
    rhoP1 = log(prob(3,:)./prob(4,:));    
    rhoP2 = log(prob(3,:)./prob(5,:));  
    rho_max=max([rhoM2(rhoM2~=inf),rhoM1(rhoM1~=inf),rhoP1(rhoP1~=inf),rhoP2(rhoP2~=inf)]);
    wetConst=max(rho_max*10^4,10^8);  %set wet cost        
    rhoP1(isnan(rhoP1)) = wetConst;    
    rhoM1(isnan(rhoM1)) = wetConst;
    rhoP2(isnan(rhoP2)) = wetConst;    
    rhoM2(isnan(rhoM2)) = wetConst;
    rhoP1(rhoP1==inf) = wetConst;    
    rhoM1(rhoM1==inf) = wetConst;
    rhoP2(rhoP2==inf) = wetConst;    
    rhoM2(rhoM2==inf) = wetConst;  
    %embed
    cover_1D=reshape(cover,N,1);
    rhoP1_1D=reshape(rhoP1,N,1);
    rhoM1_1D=reshape(rhoM1,N,1);
    rhoP2_1D=reshape(rhoP2,N,1);
    rhoM2_1D=reshape(rhoM2,N,1);
    rho0=zeros(size(rhoP1_1D));
    cost=[rhoP1_1D,rhoM1_1D,rho0,rhoP2_1D,rhoM2_1D];
    cost=cost';    
    if strcmp(coding_method,'Sim')==1   %simulator
        [lambda,pp]=calc_lambda_curve_forward_lt1(cost,mlen,N);
        pChangeP1 = pp(1,:)';
        pChangeM1 = pp(2,:)';
        pChangeP2 = pp(4,:)';
        pChangeM2 = pp(5,:)';
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
        randChange = rand(size(cover_1D));
        stego_1D = cover_1D;
        stego_1D(randChange < pChangeP2) = stego_1D(randChange < pChangeP2) + 2;
        stego_1D(randChange >= pChangeP2 & randChange < pChangeP2+pChangeP1) = stego_1D(randChange >= pChangeP2 & randChange < pChangeP2+pChangeP1) + 1;
        stego_1D(randChange >= pChangeP2+pChangeP1 & randChange < pChangeP2+pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP2+pChangeP1 & randChange < pChangeP2+pChangeP1+pChangeM1) - 1;
        stego_1D(randChange >= pChangeP2+pChangeP1+pChangeM1 & randChange < pChangeP2+pChangeP1+pChangeM1+pChangeM2) = stego_1D(randChange >= pChangeP2+pChangeP1+pChangeM1 & randChange < pChangeP2+pChangeP1+pChangeM1+pChangeM2) - 2;
        stego=reshape(stego_1D,size(cover));
    elseif strcmp(coding_method,'STC')==1
        cover_stc=int32(cover_1D');
        cost_stc=single([rhoM2_1D';rhoM1_1D';rho0';rhoP1_1D';rhoP2_1D']);
        msg=randi([0,1],mlen,1);
        msg_stc=uint8(msg');
        try
            [dist, stego_stc, n_msg_bits, l] = stc_pm2_pls_embed(cover_stc, cost_stc, msg_stc, para);
            stego_1D=double(stego_stc);
            stego=reshape(stego_1D,size(cover));
%             extr_msg = stc_ml_extract(int32(stego_stc), n_msg_bits, para);
%             success=0;
%             if all(extr_msg==msg_stc(1:numel(extr_msg)))
%                 success=1;
% %                 if numel(extr_msg)~=numel(msg_stc)
% %                     success=(numel(msg_stc)-numel(extr_msg))*100;
% %                 end
%             end
        catch
%             disp(cover_img);
%             delete(cover_img);
%             return;
            mark_success=0;
            [lambda,pp]=calc_lambda_curve_forward_lt1(cost,mlen,N);
            pChangeP1 = pp(1,:)';
            pChangeM1 = pp(2,:)';
            pChangeP2 = pp(4,:)';
            pChangeM2 = pp(5,:)';
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
            randChange = rand(size(cover_1D));
            stego_1D = cover_1D;
            stego_1D(randChange < pChangeP2) = stego_1D(randChange < pChangeP2) + 2;
            stego_1D(randChange >= pChangeP2 & randChange < pChangeP2+pChangeP1) = stego_1D(randChange >= pChangeP2 & randChange < pChangeP2+pChangeP1) + 1;
            stego_1D(randChange >= pChangeP2+pChangeP1 & randChange < pChangeP2+pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP2+pChangeP1 & randChange < pChangeP2+pChangeP1+pChangeM1) - 1;
            stego_1D(randChange >= pChangeP2+pChangeP1+pChangeM1 & randChange < pChangeP2+pChangeP1+pChangeM1+pChangeM2) = stego_1D(randChange >= pChangeP2+pChangeP1+pChangeM1 & randChange < pChangeP2+pChangeP1+pChangeM1+pChangeM2) - 2;
            stego=reshape(stego_1D,size(cover));
        end
    elseif strcmp(coding_method,'SPC')==1
        cover_stc=int32(cover_1D');
        cost_stc=single([rhoM2_1D';rhoM1_1D';rho0';rhoP1_1D';rhoP2_1D']);
        msg=randi([0,1],mlen,1);
        msg_stc=uint8(msg');
        try
            [dist, stego_stc, n_msg_bits] = spc_ml_pls_embed(cover_stc, cost_stc, msg_stc, para);
            stego_1D=double(stego_stc);
            stego=reshape(stego_1D,size(cover));
        catch
            mark_success=0;
            [lambda,pp]=calc_lambda_curve_forward_lt1(cost,mlen,N);
            pChangeP1 = pp(1,:)';
            pChangeM1 = pp(2,:)';
            pChangeP2 = pp(4,:)';
            pChangeM2 = pp(5,:)';
            RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
            randChange = rand(size(cover_1D));
            stego_1D = cover_1D;
            stego_1D(randChange < pChangeP2) = stego_1D(randChange < pChangeP2) + 2;
            stego_1D(randChange >= pChangeP2 & randChange < pChangeP2+pChangeP1) = stego_1D(randChange >= pChangeP2 & randChange < pChangeP2+pChangeP1) + 1;
            stego_1D(randChange >= pChangeP2+pChangeP1 & randChange < pChangeP2+pChangeP1+pChangeM1) = stego_1D(randChange >= pChangeP2+pChangeP1 & randChange < pChangeP2+pChangeP1+pChangeM1) - 1;
            stego_1D(randChange >= pChangeP2+pChangeP1+pChangeM1 & randChange < pChangeP2+pChangeP1+pChangeM1+pChangeM2) = stego_1D(randChange >= pChangeP2+pChangeP1+pChangeM1 & randChange < pChangeP2+pChangeP1+pChangeM1+pChangeM2) - 2;
            stego=reshape(stego_1D,size(cover));            
        end
    else 
        error('Without such coding method!');
    end  
end
S_STRUCT = img;
S_STRUCT.coef_arrays{1}=stego;
jpeg_write(S_STRUCT,stego_img);


end
