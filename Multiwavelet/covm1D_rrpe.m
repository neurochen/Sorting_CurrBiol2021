function covm=covm1D_rrpe(n,flt,maxlevel)
%covm=covm1D_rrpe(n,flt,maxlevel)
%
%  This function computes covariance matrix of 1-dimensional multiwavelet 
%  transform with oversampled ("repeated row") preprocessing. Boundaries are 
%  handled by periodic extension. For details see [SW].
%
%  Input:                                                      
%    maxlevel    integer, number of levels of decomposition    
%    flt         string of characters, name of the filter;
%                for possible names and short descriptions see coef.m  
%    n           integer, length of the intended input signal 
%                n must be of the form integer*2^maxlevel;
%                r is the number of scaling functions  
%
%  Output:                      
%    covm        r*n by r*n real array, covariance matrix of the transform
%  
%  Example of Usage:
%   covm=covm1D_rrpe(256,'ghm',3)

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);
nf=length(L(:,1));

vs=diag(ones(1,n));

covm=zeros(nf*n,n);
for i=1:n
  fp=prep1D_rr(vs(i,:),flt);
  transf=dec1D_pe(fp,flt,maxlevel);
  for j=1:n
    covm(nf*(j-1)+1:nf*j,i)=transf(:,j);
  end
end

covm=covm*covm';
