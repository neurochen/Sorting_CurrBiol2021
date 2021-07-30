function transfth=thresh_decor1D(transf,maxlevel,thresh,covname,threshtype)
%transfth=thresh_decor1D(transf,maxlevel,thresh,covname,threshtype)
%
%  This subroutine performs scalar thresholding with decorrelation of 
%  1-dimensional wavelet transform. Scaling coefficients stay untouched. 
%  Wavelet coefficients are decorrelated using diagonal blocks of the 
%  covariance matrix of the transform. After that all obtained coefficients 
%  are thresholded and correlated again.                            
%
%  Input:                                                       
%    threshtype  string of characters, either 'hard' or 'soft';
%                  threshtype='hard' invokes hard thresholding,
%                  threshtype='soft' invokes soft thresholding;
%                for details see [DJ]
%    covname     string of characters, name of covariance matrix used for 
%                decorrelation; for possible names and short descriptions see
%                coef_dcov1D.m     
%    thresh      real, value of the threshold
%    maxlevel    integer, number of levels of wavelet decomposition in transf,
%                maxlevel < 7
%    transf      r by n real array, wavelet transform of a signal;
%                r is the number of scaling functions,
%                n must be of the form integer*2^maxlevel,
%                for structure of transf see [SW]
%
%  Output:
%    transfth    r by n real array, thresholded wavelet transform
%
%  Example of Usage:
%    transfth=thresh_decor1D(transf,5,0.75,'ghmap','hard')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

n0=length(transf(1,:));
D=coef_dcov1D(covname);
r=length(D(:,1));
transfth=transf;


for i=1:maxlevel,
  d=D(:,(i-1)*r+1:i*r);
  dd=inv(sqrtm(d));
  transfth(:,n0-(2^i-1)*n0/2^i+1:n0-(2^(i-1)-1)*n0/2^(i-1))=dd*transfth(:,n0-(2^i-1)*n0/2^i+1:n0-(2^(i-1)-1)*n0/2^(i-1));
end

transfth=thresh(transfth,maxlevel,thresh,threshtype);

for i=1:maxlevel,
  d=D(:,(i-1)*r+1:i*r);
  dd=sqrtm(d);
  transfth(:,n0-(2^i-1)*n0/2^i+1:n0-(2^(i-1)-1)*n0/2^(i-1))=dd*transfth(:,n0-(2^i-1)*n0/2^i+1:n0-(2^(i-1)-1)*n0/2^(i-1));
end



