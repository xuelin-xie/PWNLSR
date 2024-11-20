# Remark
<font color="red">***Codes will be uploaded after the paper is accepted.***</font>

# PWNLSR
This is the matlab code of paper "Pixel-Weighted Nonlocal Low-Rank Subspace Representation for Hyperspectral Image Denoising".

# Introduction
Our motivation is to distinguish high-quality and low-quality pixels through weight tensors.  
![image](https://github.com/xuelin-xie/PWNLSR/blob/main/image/Flowchart.png)

# Contents 
1. Pixel-Weighted HSI (PWHSI) denoising framework:  
![image](https://github.com/xuelin-xie/PWNLSR/blob/main/image/eq1.png)
2. PWNLSR model:  
![image](https://github.com/xuelin-xie/PWNLSR/blob/main/image/eq2.png)

# Why does the PWHSI denoising framework work?
Pixel weighting mainly affects the iteration variables $\mathcal{T}$.  
![image](https://github.com/xuelin-xie/PWNLSR/blob/main/image/eq3.png)  
Without pixel weighting, all weight values default to 1, failing to dynamically and effectively balance the contributions of the noisy and restored images. By introducing pixel weighting, the auxiliary variables are more accurately estimated, thereby enhancing the performance of subsequent subspace decomposition and non-local mean operations.

# Experiments
*** Comparison Methods  
%   1.  LRMR,                     2014  TGRS  
%   2.  LLRT,                     2017  CVPR, Long time  
%   3.  FastHyDe,                 2018  J-STARS  
%   4.  LRTDGS,                   2020  TCybs  
%   5.  E-3DTV,                   2020  TIP  
%   6.  FGSLR,                    2022  TGRS, Long time  
%   7.  NGmeet,                   2022  TPAMI  
%   8.  RCTV,                     2022  TGRS  
%   9.  NS3R,                     2023  J-STARS  
%   10. TPTV,                     2023  TGRS  
%   11. 3DCTV-RPCA,               2023  TPAMI  
%   12. Ours (PWNLSR)  

API of all methods are list in "demo_PWNLSR.m"   
Run   "demo_PWNLSR.m"  to test the code for all simulation experiments.  
Run   "Real_demo.m"   to test the code with all real data.  

Theoretically, the pixel-weighted HSI (PWHSI) denoising framework can handle various types of noise, and may be more effective for impulse noise (See our code for more explorations). However, limited by the combination of NL and SR, the PWNLSR method is not particularly good at handling impulse noise or heavy stripe noise.  

# How to cite
If you use the code, please cite our papers: "Pixel-Weighted Nonlocal Low-Rank Subspace Representation for Hyperspectral Image Denoising. ..."

# Contact
If you find any bug or have any questions, please contact Xuelin Xie ( xl.xie@whu.edu.cn).
