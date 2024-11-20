clear;
clc;
%%%% Important: For standardization, we convert the original data values ​​to the range of 0-1.
% addpath(genpath(pwd));

%% Real Data
% % test1        Indian_pines   145*145*220
% load indian_pines   
% O_Img = indian_pines;

% % test2        HYDICE Urban 307*307*210
load Urban

%% Comparison Methods
%   1.  LRMR,            2014  TGRS
%   2.  LLRT,            2017  CVPR, Long time
%   3.  FastHyDe,        2018  J-STARS
%   4.  LRTDGS,          2020  TCybs
%   5.  E-3DTV,          2020  TIP
%   6.  FGSLR,           2022  TGRS
%   7.  NGmeet,          2022  TPAMI
%   8.  RCTV,            2022  TGRS
%   9.  NS3R,            2023  J-STARS
%   10. TPTV,            2023  TGRS
%   11. 3DCTV-RPCA,      2023  TPAMI, Jiangjun Peng,  Indian_pines best but long time
%   12. Ours (PWNLSR)

Comparison_Nosiy     = 1;
Comparison_LRMR      = 0;
Comparison_LLRT      = 0;   % Long time
Comparison_FastHyDe  = 0;
Comparison_LRTDGS    = 0;
Comparison_E3DTV     = 0;
Comparison_FGSLR     = 0;   % Long time
Comparison_NGmeet    = 0;
Comparison_RCTV      = 0;
Comparison_NS3R      = 0;
Comparison_TPTV      = 0;   % Long time
Comparison_3DCTV     = 0;
Comparison_PWNLSR    = 1;

%% Basic settings
i=0;
setdemorandstream(pi);  % Setting the random seed
[M, N, bands] = size(O_Img);
% Convert to a value between 0 and 1
O_Img_normalized = zeros(M, N, bands);
for b = 1:bands
    min_val = min(min(O_Img(:,:,b)));  % Minimum value of the current band
    max_val = max(max(O_Img(:,:,b)));  % Maximum value of the current band
    O_Img_normalized(:,:,b) = (O_Img(:,:,b) - min_val) / (max_val - min_val);
end
O_Img = O_Img_normalized;

%% Noise selection
nSig = 0.01/255;  % Default Value
Nosiy_case='Real data';
N_Img = O_Img;

%% Noisy
if Comparison_Nosiy == 1
    i=i+1;
    Time(i) = 0;
    Methods{i} = 'Noisy';
end

%% LRMR, 2014 TGRS
if Comparison_LRMR == 1
    i=i+1;
    rank = 4;       % Indian_pines 4,     Urban 4
    Mvalue = 20;    % Indian_pines 20,    Urban  20
    s = 4000;      % Indian_pines 4000,   Urban  4000
    stepszie = 4;
    tic
    LRMR_Img = LRMR_HSI_denoise(N_Img, rank, Mvalue, s, stepszie);
    Time(i) = toc;
    Methods{i} = 'LRMR';
end

%% LLRT, 2017 CVPR, Long time
if Comparison_LLRT == 1
    i=i+1;
    Par_LLRT = ParSetC(nSig*255, bands);
    tic
    LLRT_Img = LLRT_DeNoising(N_Img*255, O_Img, Par_LLRT)/255;
    Time(i) = toc;
    Methods{i} = 'LLRT';
end

%% FastHyDe, 2018 J-STARS
if Comparison_FastHyDe == 1
    i=i+1;
    noise_type = 'additive';
    iid = 1;
    k_subspace = 6;       % Indian_pines 4,  Urban 6
    tic
    FastHyDe_Img = FastHyDe(N_Img, noise_type, iid, k_subspace); % since there exists timing perform in FastHyDe
    Time(i) = toc;
    Methods{i} = 'FastHyDe';
end

%% LRTDGS, 2020 TCybs
if Comparison_LRTDGS == 1
    i=i+1;
    lambda1 = 0.8;               % parameter for group sparse term  Indian_pines 0.8,
    lambda2 = 50/(sqrt(M*N));     % parameter for sparse noise
    Rank = [120,120,10];          % parameter for Tucker rank
    tic
    [LRTDGS_Img,~,~,~] = LRTDGS_Denoising(N_Img,lambda1,lambda2,Rank);
    Time(i) = toc;
    Methods{i} = 'LRTDGS';
end

%% E-3DTV, 2020 TIP
if Comparison_E3DTV == 1
    i=i+1;
    rank = [13, 13, 13];    % PaC [13, 13, 13], CAVE [13, 13, 8], WDC [13, 8, 8]
    tau = 0.004 * sqrt(M*N);   % Indian_pines 0.04 * sqrt(M*N) WDC 0.003 * sqrt(M*N);
    tic
    E3DTV_Img = EnhancedTV(N_Img, tau, rank);  % since there exists timing perform in E3DTV
    Time(i) = toc;
    Methods{i} = 'E-3DTV';
end

%% FGSLR, 2022 TGRS
if Comparison_FGSLR == 1
    i=i+1;
    opt.alpha = 1;
    opt.beta = 0.5;          % TV regularization:set as 0.1 or 0.5            CAVE 0.1  WDC
    opt.lambda = 0.01;       % sparse noise
    opt.mu = 5;              % penalty parameter:set as 5 or 10      PaC: 5   CAVE: 15  WDC
    opt.rho = 0.1;           % proximal parameter
    opt.delta = 5;           % noise parameter:set as 0.5 or 5                CAVE: 15  WDC
    opt.r = 20;
    opt.regul_B = 'L21';
    tic
    [X, A, B] = FGSLR_TV_PAM(N_Img, O_Img, opt);
    FGSLR_Img = reshape(X', M, N, bands);
    Time(i) = toc;
    Methods{i} = 'FGSLR';
end

%% NGmeet, 2022 TPAMI     also 2019 CVPR
if Comparison_NGmeet == 1
    i=i+1;
    noiselevel = nSig*ones(1,80);
    ParNG = ParSetH(255*mean(noiselevel),bands);
    tic
    NGmeet_Img = NGmeet_DeNoising(255*N_Img, ParNG)/255;  %NGmeet denoisng function
    Time(i) = toc;
    Methods{i} = 'NGmeet';
end

%% RCTV, 2022 TGRS
if Comparison_RCTV == 1
    i=i+1;
    r=5;         % CAVE 13,  PaC 5,  WDC 5
    beta = 10;
    lambda = 5;  % 5,0.5   % Indian_pines 1,  PaC 1.5,  CAVE 0.6,   WDC  1.5
    tau = [0.7, 0.95];  % WDC [0.7, 0.85]  Indian_pines [0.7,0.95]    PaC CAVE [0.99, 0.85]
    tic
    RCTV_Img = RCTV(N_Img, beta, lambda, tau, r);
    Time(i) = toc;
    Methods{i} = 'RCTV';
end

%% NS3R, 2023 J-STARS
if Comparison_NS3R == 1
    i=i+1;
    Par_NS3R = ParSetNS3R(nSig, size(N_Img));
    tic
    NS3R_Img = NS3R(N_Img, Par_NS3R);
    Time(i) = toc;
    Methods{i} = 'NS3R';
end

%% TPTV, 2023 TGRS
if Comparison_TPTV == 1
    i=i+1;
    param.Rank = [7,7,5];           % Indian_pines, PaC: [7,7,5],  CAVE [10,10,9]  WDC [7,7,5]
    param.initial_rank = 2;
    param.maxIter = 50;
    param.lambda  = 4e-3*sqrt(M*N);  % 0.1;   % WDC, Indian_pines: 4e-3*sqrt(M*N),  PaC: 3e-3*sqrt(M*N),  cave: 6e-3*sqrt(M*N)
    tic
    [output_image,U_x,V_x,E] = WETV(N_Img, O_Img, param);
    TPTV_Img = reshape(output_image, M, N, bands);
    Time(i) = toc;
    Methods{i} = 'TPTV';
end

%% 3DCTV_rpca, 2023 TPAMI
if Comparison_3DCTV == 1
    i=i+1;
    tic
    weight=0.5;
    CTV3D_Img = ctv_rpca(N_Img, weight);
    Time(i) = toc;
    Methods{i} = '3DCTV';
end
% imshow(mat2gray(cat(3,CTV3D_Img(:,:,105),CTV3D_Img(:,:,160),CTV3D_Img(:,:,220))));

%% PWNLSR, Ours
if Comparison_PWNLSR == 1
    i=i+1;
    ParPWNLSR_real=ParSetPWNLSR_real(nSig,bands);
    tic
    [PWNLSR_Img, W] = PWNLSR(N_Img, ParPWNLSR_real);
    Time(i) = toc;
    Methods{i} = 'Ours (PWNLSR)';
end

% Indian_pines
% figure; imshow(mat2gray(cat(3,N_Img(:,:,105),N_Img(:,:,160),N_Img(:,:,220))));
% figure; imshow(mat2gray(cat(3,PWNLSR_Img(:,:,105),PWNLSR_Img(:,:,160),PWNLSR_Img(:,:,220))));

% Urban
% figure; imshow(mat2gray(cat(3,N_Img(:,:,3),N_Img(:,:,207),N_Img(:,:,119))));
% figure; imshow(mat2gray(cat(3,PWNLSR_Img(:,:,3),PWNLSR_Img(:,:,207),PWNLSR_Img(:,:,119))));