function thadj=coef_thadj(pflt)
%thadj=coef_thadj(pflt)
%
%  This function returns matrix of mean variances of a given transform.
%  Each element of this matrix corresponds to certain length of the input
%  signal and certain depth of wavelet decomposition. For details see [SW].
%
%  Input:                                                   
%    pflt        string of characters, name of the transform; for admissible
%                names and short descriptions see below.
%
%  Output: 
%    thadj        7 by 7 real array, mean variances of the transform;
%                 data in thadj is organized as follows:
%                   each row of thadj corresponds to the certain length
%                   of the signal (32, 64, ..., 2048, from top to bottom),
%                   each column of thadj corresponds to certain number of
%                   levels of wavelet decomposition (2, 3, ..., 7, from
%                   left to right)
%
%  Admissible Names of the Transforms:                              
%  (all names of wavelet transforms are from coef.m,
%   all names of prefilters are from coef_prep.m)
% 
%    'ghmrr'      'ghm' multiwavelet transform with oversampled preprocessing,
%    'ghmghmap'   'ghm' multiwavelet transform with with 'ghmap' prefiltering,
%    'ghmghmorap' 'ghm' multiwavelet transform with with 'ghmorap' 
%                 prefiltering,
%    'clrr'       'cl' multiwavelet transform with oversampled preprocessing,
%    'clclap'     'cl' multiwavelet transform with with 'clap' prefiltering,
%    'sa4rr'      'sa4' multiwavelet transform with oversampled preprocessing,
%    'sa4sa4ap'   'sa4' multiwavelet transform with 'sa4ap' prefiltering,
%    'd4scalar'   'd4' scalar wavelet transform,
%    'la8scalar'  'la8' scalar wavelet transform,
%    'bih54nrr'   'bih54n' multiwavelet transform with oversampled 
%                 preprocessing,
%    'bih54nbih5ap'  'bih54n' multiwavelet transform with 'bih5ap' prefiltering,
%    'bih52srr'   'bih52s' multiwavelet transform with oversampled 
%                 preprocessing,
%    'bih52sbih5ap'  'bih54n' multiwavelet transform with 'bih5ap' prefiltering,
%    'bi9scalar' 'bi9' scalar biorthogonal wavelet transform, 
% 
%
%  Example of Usage:
%   thadj=coef_thadj('bih552sbih5ap')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

a=ones(7);

if strcmp(pflt,'ghmrr')
  thadj=sqrt(3/4)*a;
elseif strcmp(pflt,'ghmghmap')
  thadj=0.980*a;
elseif strcmp(pflt,'ghmghmorap')
  thadj=a;
elseif strcmp(pflt,'clrr')
  thadj=sqrt(1/2)*a;
elseif strcmp(pflt,'clclap')
  thadj=0.3711*a;
elseif strcmp(pflt,'sa4rr')
  thadj=sqrt(1/2)*a;
elseif strcmp(pflt,'sa4sa4ap')
  thadj=a;
elseif strcmp(pflt,'d4scalar')
  thadj=a;
elseif strcmp(pflt,'la8scalar')
  thadj=a;
elseif strcmp(pflt,'bih54nrr')
  thadj=[2.4993    0         0         0         0         0         0;
         2.4685    3.0360    0         0         0         0         0;
         2.5113    3.1514    3.4429    0         0         0         0;
         2.4935    3.1412    3.5917    4.0982    0         0         0;
         2.4887    3.0965    3.6248    4.1330    4.4911    0         0;
         2.4946    3.0963    3.6439    4.0049    4.5797    4.7930    0;
         2.4867    3.0894    3.6174    4.0308    4.5419    5.0452    5.3843];
elseif strcmp(pflt,'bih54nbih5ap')
  thadj=[1.4128    0         0         0         0         0         0;
         1.4379    1.7458    0         0         0         0         0;
         1.4624    1.7665    2.0159    0         0         0         0;
         1.4913    1.8085    2.0560    2.2913    0         0         0;
         1.4630    1.8159    2.0879    2.3826    2.5088    0         0;
         1.4885    1.8060    2.0923    2.3942    2.5418    2.8059    0;
         1.4873    1.8235    2.0989    2.3659    2.6201    2.9911    3.1561];
elseif strcmp(pflt,'bih52srr')
  thadj=[3.0722    0         0         0         0         0         0;
         3.0190    2.7322    0         0         0         0         0;
         3.0945    2.8258    2.3783    0         0         0         0;
         3.0807    2.8146    2.4429    2.0525    0         0         0;
         3.0708    2.8130    2.3904    2.1213    1.9204    0         0;
         3.0717    2.7951    2.4074    2.1107    1.9203    1.8158    0;
         3.0648    2.8064    2.3992    2.0959    1.9120    1.8181    1.7606];
elseif strcmp(pflt,'bih52sbih5ap')
  thadj=[1.4009    0         0         0         0         0         0;
         1.3620    1.2982    0         0         0         0         0;
         1.4346    1.3634    1.1865    0         0         0         0;
         1.4471    1.3677    1.1765    1.0112    0         0         0;
         1.4244    1.3742    1.1962    1.0403    0.9428    0         0;
         1.4441    1.3569    1.1864    1.0420    0.9469    0.8991    0;
         1.4474    1.3706    1.1872    1.0386    0.9473    0.8975    0.8655];
elseif strcmp(pflt,'bi9scalar')
  thadj=[1.0130    0         0         0         0         0         0;
         1.0055    1.0136    0         0         0         0         0;
         1.0087    1.0123    1.0113    0         0         0         0;
         1.0093    1.0114    1.0109    1.0127    0         0         0;
         1.0096    1.0113    1.0112    1.0120    1.0126    0         0;
         1.0098    1.0119    1.0120    1.0130    1.0129    1.0129    0;
         1.0095    1.0118    1.0126    1.0123    1.0129    1.0131    1.0128];
end
