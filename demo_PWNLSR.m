clear;
clc;
%%%% Important: For standardization, we convert the original data values ​​to the range of 0-1.
% addpath(genpath(pwd));

%% Test Data
% % test1 PaC, University of Pavia, Cropped image 300*300*103
load paviaU_cropped
O_Img = paviaU_cropped;

% % test2,  WDC, Washington DC MALL, 256*256*191 
% load Ori_WDC
% O_Img = Img;

% % test 3,  Salinas, 217*217*224
% load salinas

%% Comparison Methods
%   1.  LRMR,            2014  TGRS
%   2.  LLRT,            2017  CVPR, Long time
%   3.  FastHyDe,        2018  J-STARS
%   4.  LRTDGS,          2020  TCybs
%   5.  E-3DTV,          2020  TIP
%   6.  FGSLR,           2022  TGRS, Long time
%   7.  NGmeet,          2022  TPAMI
%   8.  RCTV,            2022  TGRS
%   9.  NS3R,            2023  J-STARS
%   10. TPTV,            2023  TGRS
%   11. 3DCTV-RPCA,      2023  TPAMI, Jiangjun Peng
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
Nosiy_case=5;
if Nosiy_case==1
    % Case1.  Gaussian noise  0.0392
    nSig = 10/255;  % Noise level
    N_Img_Gaussian = O_Img + nSig * randn(size(O_Img));
    N_Img = N_Img_Gaussian;
elseif Nosiy_case==2
    % Case2.  Gaussian noise  0.1176
    nSig = 30/255;  % Noise level
    N_Img_Gaussian = O_Img + nSig * randn(size(O_Img));
    N_Img = N_Img_Gaussian;
elseif Nosiy_case==3
   % Case3.  Gaussian noise  0.1961
    nSig = 50/255;  % Noise level
    N_Img_Gaussian = O_Img + nSig * randn(size(O_Img));
    N_Img = N_Img_Gaussian;
elseif Nosiy_case==4
    % Case 4. Gaussian noise  0.3922
    nSig = 100/255;  % Noise level
    N_Img_Gaussian = O_Img + nSig * randn(size(O_Img));
    N_Img = N_Img_Gaussian;
elseif Nosiy_case==5
    % Case 5. Mixed noise 1: Gaussian 0.5 + stripe noise [0.2, 0.1]
    nSig = 50/255;
    N_Img_Gaussian = O_Img + nSig * randn(size(O_Img));
    N_Img_Mix1 = N_Img_Gaussian;
    stripeRatio = 0.2;  % Stripe ratio
    numStripes = round(N * stripeRatio);  % Number of stripes
    stripeOffset = 0.1;  % Fixed offset of stripes
    for b = 1:bands
        stripeColumns = randperm(N, numStripes);
        for c = stripeColumns
            N_Img_Mix1(:, c, b) = N_Img_Mix1(:, c, b) + stripeOffset;
        end
    end
    N_Img =N_Img_Mix1;
end

%% Noisy
if Comparison_Nosiy == 1
    i=i+1;
    Time(i) = 0;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, N_Img*255);
    disp(['Method Name: None', ', Time = None', ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'Noisy';
end

%% LRMR, 2014 TGRS
if Comparison_LRMR == 1
    i=i+1;
    rank = 5;                     
    Mvalue = 25;                           
    s = 4000;                      
    stepszie = 4;
    tic
    LRMR_Img = LRMR_HSI_denoise(N_Img, rank, Mvalue, s, stepszie);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, LRMR_Img*255);
    disp(['Method Name: LRMR', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'LRMR';
end

%% LLRT, 2017 CVPR, Long time
if Comparison_LLRT == 1
    i=i+1;
    Par_LLRT = ParSetC(nSig*255, bands);
    tic
    LLRT_Img = LLRT_DeNoising(N_Img*255, O_Img, Par_LLRT)/255;
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, LLRT_Img*255);
    disp(['Method Name: LLRT', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'LLRT';
end

%% FastHyDe, 2018 J-STARS
if Comparison_FastHyDe == 1
    i=i+1;
    noise_type = 'additive';
    iid = 1;
    k_subspace = 6;       
    tic
    FastHyDe_Img = FastHyDe(N_Img, noise_type, iid, k_subspace); % since there exists timing perform in FastHyDe
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, FastHyDe_Img*255);
    disp(['Method Name: FastHyDe', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'FastHyDe';
end

%% LRTDGS, 2020 Tcybs
if Comparison_LRTDGS == 1
    i=i+1;
    lambda1 = 0.052;                              %  PaC, Salinas: 0.052;    WDC: 0.05; 
    lambda2 = 50/(sqrt(M*N));   
    Rank_LRTDGS = [140,140,4];                    %  PaC: [140,140,4];  WDC:[180,180,4];  Salinas: [200,200,5]
    tic
    LRTDGS_Img = LRTDGS_Denoising(N_Img, lambda1,lambda2,Rank_LRTDGS);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, LRTDGS_Img*255);
    disp(['Method Name: LRTDGS', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'LRTDGS';
end

%% E-3DTV, 2020 TIP
if Comparison_E3DTV == 1
    i=i+1;
    rank = [13, 8, 8];
    tau = 0.0019 * sqrt(M*N);
    tic
    E3DTV_Img = EnhancedTV(N_Img, tau, rank); 
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, E3DTV_Img*255);
    disp(['Method Name: E3DTV', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'E-3DTV';
end

%% FGSLR, 2022 TGRS
if Comparison_FGSLR == 1
    i=i+1;
    opt.alpha = 1;
    opt.beta = 0.5;                           % PaC, WDC: 0.5;     Salinas: 0.1
    opt.lambda = 0.01;         
    opt.mu = 10;               
    opt.rho = 0.1;         
    opt.delta = 0.5;           
    opt.r = 20;
    opt.regul_B = 'L21';
    tic
    [X, A, B] = FGSLR_TV_PAM(N_Img, O_Img, opt);
    FGSLR_Img = reshape(X', M, N, bands);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, FGSLR_Img*255);
    disp(['Method Name: FGSLR', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'FGSLR';
end

%% NGmeet, 2022 TPAMI     also 2019 CVPR
if Comparison_NGmeet == 1
    i=i+1;
    noiselevel = nSig*ones(1,80);
    ParNG = ParSetH(255*mean(noiselevel),bands);
    tic
    NGmeet_Img = NGmeet_DeNoising(255*N_Img, ParNG)/255; 
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, NGmeet_Img*255);
    disp(['Method Name: NGmeet', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'NGmeet';
end

%% RCTV, 2022 TGRS
if Comparison_RCTV == 1
    i=i+1;
    r=3;                   
    beta = 10;
    lambda = 1.5;    
    tau = [0.9, 0.7]; 
    tic
    RCTV_Img = RCTV(N_Img, beta, lambda, tau, r);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, RCTV_Img*255);
    disp(['Method Name: RCTV', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'RCTV';
end

%% NS3R, 2023 J-STARS
if Comparison_NS3R == 1
    i=i+1;
    Par_NS3R = ParSetNS3R(nSig, size(N_Img));
    tic
    NS3R_Img = NS3R(N_Img, Par_NS3R);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, NS3R_Img*255);
    disp(['Method Name: NS3R', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'NS3R';
end

%% TPTV, 2023 TGRS
if Comparison_TPTV == 1
    i=i+1;
    param.Rank = [7,7,5];           
    param.initial_rank = 2;
    param.maxIter = 50;
    param.lambda  = 3e-3*sqrt(M*N);             % PAC, Salinas: 3e-3*sqrt(M*N);   WDC: 4e-3*sqrt(M*N)
    tic
    [output_image,U_x,V_x,E] = WETV(N_Img, O_Img, param);
    TPTV_Img = reshape(output_image, M, N, bands);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, TPTV_Img*255);
    disp(['Method Name: TPTV', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'TPTV';
end

%% 3DCTV_pca, 2023 TPAMI
if Comparison_3DCTV == 1
    i=i+1;
    tic
    CTV3D_Img = ctv_pca(N_Img, nSig);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, CTV3D_Img*255);
    disp(['Method Name: 3DCTV', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = '3DCTV';
end

%% PWNLSR, Ours
if Comparison_PWNLSR == 1
    i=i+1;
    bands=size(N_Img,3);
    ParPWNLSR=ParSetPWNLSR(nSig,bands);
    tic
    [PWNLSR_Img, W] = PWNLSR(N_Img, ParPWNLSR);
    % W = reshape(10000*W,M*N, bands);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i), ERGAS(i), MSAM(i)] = MSIQA(O_Img*255, PWNLSR_Img*255);
    disp(['Method Name: Ours (PWNLSR)', ', Time = ' num2str(Time(i)), ', MPSNR = ' num2str(PSNR(i),'%5.2f')  ...
        ', MSSIM = ' num2str(SSIM(i),'%5.4f'), ', The case of noise is: ' num2str(Nosiy_case)]);
    Methods{i} = 'Ours (PWNLSR)';
end

%% Save the results
Indexes = [PSNR; SSIM; FSIM; ERGAS; MSAM; Time];
results = array2table(Indexes, ...
    'VariableNames', Methods', ...
    'RowNames', {'PSNR', 'SSIM', 'FSIM', 'ERGAS', 'MSAM', 'Time'});
results.Properties.DimensionNames = {'Indexes', 'Methods'};
disp(results);

save_results=0;  % save=1, Write to the excel table
if save_results==1
    writetable(results, '表1_PaC_Results(sigma=50+[0.2,0.1]).csv', 'WriteRowNames', true);
end

% % PaC image bands
% imshow(mat2gray(cat(3,N_Img(:,:,10),N_Img(:,:,30),N_Img(:,:,60))));  
% imshow(mat2gray(cat(3,PWNLSR_Img(:,:,10),PWNLSR_Img(:,:,30),PWNLSR_Img(:,:,60))));
