function D=coef_dcov1D(pflt)
%D=coef_dcov1D(pflt)
%
%  This function returns diagonal blocks of the covariance matrix 
%  corresponding to 6 levels of 1-dimensional multiwavelet transform with 
%  preprocessing. For details see [SW].//? ≤√¥“‚Àº?
%  
%  [SW]  V. Strela and A. T. Walden, "Signal and Image Denoising via Wavelet 
%        Thresholding: Orthogonal and Biorthogonal, Scalar and Multiple Wavelet 
%        Transforms", Imperial College, Statistics Section, 
%        Technical Report TR-98-01 (1998).
%
%
%  Input:                                                        
%    pflt       string of characters, name of transform used to generate 
%               the covariance matrix; for admissible names see below 
%
%  Output:                                                       
%    D          r by 6*r real array, diagonal blocks of the covariance matrix;
%               r is the number of scaling functions; 
%               blocks in D are organized as follows: D=[D1 D2 ... D6], 
%               Dk corresponds to the k-th level of the transform 
%
%  Admissible Names of the Transforms: 
%  (names of wavelet transforms are from coef.m; 
%   names of prefilters are from coef_prep.m)  
%   'ghmap'     'ghm' multiwavelet transform with 'ghmap' prefilter
%   'ghmorap'   'ghm' multiwavelet transform with 'ghmorap' prefilter 
%   'ghmrr'     'ghm' multiwavelet transform with oversampled preprocessing   
%   'clap'      'cl' multiwavelet transform with 'clap' prefilter  
%   'clrr'      'cl' multiwavelet transform with oversampled preprocessing   
%   'bih52sap'  'bih52s' multiwavelet transform with 'bih5ap' prefilter 
%   'bih52srr'  'bih52s' multiwavelet transform with oversampled preprocessing
%   'sa4ap'     'sa4' multiwavelet transform with 'sa4ap' prefilter  
%   'bi9'       'bi9' wavelet transform without preprocessing
%
%  Example of Usage:
%   D=coef_dcov1D('ghmorap')                

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

if strcmp(pflt,'ghmap')
  D=[0.5626 0      1.0540 0      1.3518 0      1.4853 0      1.5262 0      1.5377 0
     0      0.6863 0      1.2346 0      1.4192 0      1.5078 0      1.5321 0      1.5393];

elseif strcmp(pflt,'ghmorap')
  D=[ 1 0 1 0 1 0 1 0 1 0 1 0;
      0 1 0 1 0 1 0 1 0 1 0 1];

elseif strcmp(pflt,'ghmrr')
  D=[0.2700 0.2121  1.2683 -0.0212 2.5251 0.0021  2.8796 -0.0002 2.9698 0      2.9925 0
     0.2121 1.3400 -0.0212  1.9486 0.0021 2.7117 -0.0002  2.9269 0      2.9817 0      2.9954]/2;

elseif strcmp(pflt,'clap')
  D=[0.1314 -0.0000 0.1408 0      0.1281 0      0.1264 0      0.1253 0      0.1251 0
    -0.0000  0.1489 0      0.1453 0      0.1305 0      0.1271 0      0.1256 0      0.1252];

elseif strcmp(pflt,'clrr') 
  D=[0.7500 0      0.3815 0      0.8782 0      0.9439 0      0.9869 0      0.9955 0
     0      0.0625 0      0.2030 0      0.7822 0      0.9182 0      0.9780 0      0.9931];

elseif strcmp(pflt,'bih52sap')
  D=[0.0547 0      0.2030 0      0.3226 0      0.3825 0      0.4057 0      0.4136 0
     0      0.2988 0      1.5711 0      2.5777 0      3.0326 0      3.1989 0      3.2549];

elseif strcmp(pflt,'bih52srr')
  D=[0.1875 0      0.3174 0      1.8084 0       2.7501 0       3.1415 0       3.2756 0
     0      0.1406 0      0.8370 0      14.6939 0      22.0676 0      24.8666 0      25.7993];

elseif strcmp(pflt,'sa4ap')
  D=[ 1 0 1 0 1 0 1 0 1 0 1 0;
      0 1 0 1 0 1 0 1 0 1 0 1];

elseif strcmp(pflt,'bi9')
  D=[0.9830 1.1186 1.0443 1.0037 0.9901 0.9861];

end





