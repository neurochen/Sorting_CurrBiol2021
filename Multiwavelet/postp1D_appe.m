function fhat=postp1D_appe(fhatp,pflt)
%fhat=postp1D_appe(fhatp,pflt)
%
%  This function performs critically sampled postprocessing of reconstructed 
%  preprocessed 1-dimensional signal. Boundaries are handled by periodic 
%  extension. For description of the algorithm see [SW].
%
%  这是后滤波；这个函数与prep1D_appe.m共同使用；
%
%  Input:                                                    
%    pflt       string of characters, name of the prefilter;
%               for possible names and short descriptions see coef_prep.m 
%    fhatp      r by n/r real array, reconstructed preprocessed data; 
%               r is the number of scaling functions; each row corresponds 
%               to a scaling function
%
%  Output:                    
%    fhat       1 by n real array, reconstructed signal    
%
%  Example of Usage:
%   fhat=postp1D_appe(fhatp,'clap')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[PR,PO]=coef_prep(pflt);
n=length(fhatp);

[nf,np]=size(PO);
np=np/nf;

pairs=zeros(nf,n);

for i=1:np-2+1
  for j=1:np
    if i+j-np<1
      pairs(:,i)=pairs(:,i)+PO(:,nf*(j-1)+1:nf*j)*fhatp(:,n-i-j+np);
    else
      pairs(:,i)=pairs(:,i)+PO(:,nf*(j-1)+1:nf*j)*fhatp(:,i+j-np);
    end
  end
end

for i=np:n
  for j=1:np
    pairs(:,i)=pairs(:,i)+PO(:,nf*(j-1)+1:nf*j)*fhatp(:,i+j-np);
  end
end

fhat=zeros(1,nf*n);
for i=1:n
  for j=1:nf
    fhat(nf*(i-1)+j)=pairs(j,i);
  end
end