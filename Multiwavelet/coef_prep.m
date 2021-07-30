function [PR,PO]=coef_prep(pflt)
%[PR,PO]=coef_prep(pflt)
%
%  This function returns coefficients of given pre- and postfilters to be 
%  used in critical sampled preprocessing. For details see [SW].//返回的是预滤波与后滤波算子
%
%  [SW]  V. Strela and A. T. Walden, "Signal and Image Denoising via Wavelet 
%      Thresholding: Orthogonal and Biorthogonal, Scalar and Multiple Wavelet 
%      Transforms", Imperial College, Statistics Section, 
%      Technical Report TR-98-01 (1998).
%
%  Input:                                                        
%    pflt       string of characters, name of the prefilter; for admissible 
%               names and brief descriptions of filters see below 
%
%  Output: 
%    PR         预滤波
%               r by r*l real array, prefilter; 
%               r is the number of scaling functions, 
%               l is the number of coefficients in the prefilter;
%               PR is organized as follows: PR=[P1 P2 ... Pl] 
%
%    PO        后滤波 
%               r by r*m real array, postfilter;
%               m is the number of coefficients in the postfilter;
%               PO is organized as follows: PO=[R1 R2 ... Rl]
%
%  Admissible Names of the Prefilters and their Properties:       
%  (all names of the filter banks are from coef.m;   
%   r is the number of scaling functions,
%   l is the number of coefficients in the prefilter,                 
%   m is the number of coefficients in the postfilter,              
%   A is maximal approximation order preserved by the prefilter, see [SW])
%
%   'ghmap'     biorthogonal interpolation prefilter for the 'ghm' multifilter 
%               bank (see [XGHS, SHSTH]); r=2, l=2, m=2, A=2            //双正交插值预滤波
%
%  [XGHS] X.-G. Xia, J. S. Geronimo, D. P. Hardin, and B. W. Suter, 
%         "Design of prefilters for discrete multiwavelet transforms",
%          IEEE Trans. on SP, 44 (1996) 25-35.
%  [SHSTH] V. Strela, P. Heller, G. Strang, P. Topiwala, and C. Heil,
%         "The application of multiwavelet filter banks to signal and 
%          image processing", IEEE Trans. on Image Proc. (1998).
%
%   'ghmorap'   orthogonal approximating prefilter for the 'ghm' multifilter   
%               bank (see [HR]); r=2, l=2, m=2, A=2                  //正交逼近预滤波
%
%   [HR]  D. P. Hardin and D. W. Roach, "Multiwavelet prefilters I: 
%         Orthogonal prefilters preserving approximation order p <= 2",
%         preprint (1997).
%
%   'clap'      biorthogonal interpolation prefilter for 'cl' multifilter
%               bank (see [SW]); r=2, l=1, m=1, A=2
%   'bih5ap'    biorthogonal interpolation prefilter for 'bih52s' or 'bih54n'
%               multifilter bank (see [SW]); r=2, l=1, m=1, A=2
%   'bih3ap'    biorthogonal interpolation prefilter for 'bih32s'  multifilter
%               bank (see [SW]); r=2, l=1, m=1, A=2
%   'sa4ap'     orthogonal prefilter for 'sa4'  multifilter bank (see [STT]); 
%               r=2, l=1, m=1, A=1
%    [STT] L.-X. Shen, H. H. Tan, and J. Y. Tham, "Symmetric-antisymmetric 
%          orthonormal multiwavelets and related scalar wavelets", 
%          preprint (1997).
% 
%   'bighm2ap'  biorthogonal interpolation prefilter for 'bighm2' multifilter
%               bank r=2, l=1, m=1, A=2
%   'id'        orthogonal identity prefilter (is used with balanced
%               multiwavelets, see [Se]) r=2, l=1, m=1             // 正交算子预滤波，用于平衡小波
%
%    [Se]  I. Selesnick, "Cardinal Multiwavelets and the Sampling Theorem",
%          prerprint (1998).
%
%  Example of Usage:
%   [PR,PO]=coef_prep('sa4ap') 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

if strcmp(pflt,'ghmap')        %%ghn多小波用双正交插值预滤波
  PR=[3/(8*sqrt(2)) 10/(8*sqrt(2)) 3/(8*sqrt(2)) 0;
      0             0              1             0];
  PO=[0  1    0             0;
      0 -3/10 4*sqrt(2)/5 -3/10];

elseif strcmp(pflt,'ghmorap')   %%ghm多小波用正交逼近预滤波
  PR=[0.11942337067748,0.99158171438258,0.04967860804828,-0.00598315472909;
     -0.00598315472909,-0.04967860804828,0.9915817143825,-0.11942337067748];
  PO(:,1:2)=PR(:,3:4)';
  PO(:,3:4)=PR(:,1:2)';

elseif strcmp(pflt,'clap')      %%cl多小波用双正交插值预滤波
  PR=[1/4            1/4;
      1/(1+sqrt(7)) -1/(1+sqrt(7))];
  PO=[2  (1+sqrt(7))/2
      2 -(1+sqrt(7))/2];

elseif strcmp(pflt,'bih5ap')   %%双正交h5多小波用双正交插值预滤波
  PR=[1/4  1/4;
     -1 1];
  PO=[2 -1/2
      2  1/2];

elseif strcmp(pflt,'bih3ap')    %%双正交h3多小波用双正交插值预滤波
  PR=[1/4  1/4;
     -1/15 1/15];
  PO=[2 -15/2
      2  15/2];

elseif strcmp(pflt,'sa4ap')      %%sa4多小波用正交预滤波
  PR=[1,1;
      -1,1]/sqrt(2);
  PO=[1,-1;
      1,1]/sqrt(2);

elseif strcmp(pflt,'bighm2ap')    %%双正交ghm多小波用双正交插值预滤波
  PR=[1/4  1/4;
     1/3 -1/3];
  PO=[2  3/2
      2  -3/2];

elseif strcmp(pflt,'id')           %%平衡多小波，用于平衡多小波
  PR=[1 0;
      0 1];
  PO=[1 0;
      0 1];
elseif strcmp(pflt,'minial')      %% for GHM multiwavelets
    PR=[2*sqrt(2) -sqrt(2);  1           0 ];
    PO=[0           1     ;  -1/sqrt(2)  2 ];
    
elseif strcmp(pflt,'orap')       %% for GHM multiwavelets  the result seems bad 
    PR=(sqrt(6)/6)*[1-sqrt(2) 1+sqrt(2) ; 1+sqrt(2) -1+sqrt(2)];
    PO=inv(PR);
    
end
