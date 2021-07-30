function covm=covm1D_appe(n,flt,pflt,maxlevel)
%covm=covm1D_appe(n,flt,pflt,maxlevel)
%
%  This function computes covariance matrix of 1-dimensional multiwavelet 
%  transform with critically sampled preprocessing. Boundaries are handled by 
%  periodic extension. For details see [SW].
%
%  Input:                                                        
%    maxlevel    integer, number of levels of decomposition      
%    pflt        string of characters, name of the prefilter;
%                for possible names and short descriptions see coef_prep.m 
%    flt         string of characters, name of the filter;
%                for possible names and short descriptions see coef.m 
%    n           integer, length of the intended input signal;
%                n must be of the form integer*r*2^maxlevel;     
%                r is the number of scaling functions
%                                                             
%  Output:                                                    
%    covm        n by n real array, covariance matrix of the transform
%  
%  Example of Usage:
%   covm=covm1D_appe(256,'ghm','ghmap',3) 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);
nf=length(L(:,1));

vs=diag(ones(1,n));

covm=zeros(n,n);
for i=1:n
  fp=prep1D_appe(vs(i,:),pflt);
  transf=dec1D_pe(fp,flt,maxlevel);
  for j=1:round(n/nf)
    covm(nf*(j-1)+1:nf*j,i)=transf(:,j);
  end
end

covm=covm*covm';



