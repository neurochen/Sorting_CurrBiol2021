function transfth=thresh_vec1D(transf,maxlevel,thresh,covname,threshtype)
%transfth=thresh_vec1D(transf,maxlevel,thresh,covname,threshtype)
%
%  This function performs vector thresholding with decorrelation of 
%  1-dimensional multiwavelet transform. Scaling coefficients stay untouched.
%  For description of the algorithm see  [DS], [SW]
%
%  Input:   
%    threshtype  string of characters, either 'hard' or 'soft';
%                  threshtype='hard' invokes hard thresholding,
%                  threshtype='soft' invokes soft thresholding;
%                for details see [DS]
%    covname     string of characters, name of covariance matrix used for 
%                decorrelation; for possible names and short descriptions 
%                see coef_dcov1D.m
%    thresh      real, value of the threshold 
%    maxlevel    integer, number of levels of wavelet decomposition in transf;
%                maxlevel < 7 
%    transf      r by n real array, wavelet transform of a signal;
%                r is the number of scaling functions,   
%                n must be of the form integer*2^maxlevel;
%                for structure of transf see [SW]
%
%  Output:  
%    transfth    r by n real array, thresholded multiwavelet coefficients
%
%  Example of Usage:
%    transfth=thresh_vec1D(transf,5,0.75,'clap','hard')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

n0=length(transf(1,:));
D=coef_dcov1D(covname);
r=length(D(:,1));
transfth=transf;

for i=1:maxlevel,
  d=D(:,(i-1)*r+1:i*r);
  dd=inv(d);
  for j=n0-(2^i-1)*n0/2^i+1:n0-(2^(i-1)-1)*n0/2^(i-1),
    tj=sqrt(transfth(:,j)'*dd*transfth(:,j));
    if tj<thresh
      transfth(:,j)=zeros(1,r)';
    end
    if strcmp(threshtype,'soft'),
      transfth(:,j)=transfth(:,j)*(tj-thresh)/tj;
    end
  end
end






