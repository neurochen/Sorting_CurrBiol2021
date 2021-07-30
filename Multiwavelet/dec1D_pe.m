function transf=dec1D_pe(fp,flt,maxlevel)
%transf=dec1D_pe(fp,flt,maxlevel)
%
%  This function computes discrete multiwavelet transform of a preprocessed 
%  1-dimensional signal. Boundaries are handled by periodic extension.
%  For description of the algorithm see [SW].
%
%  对已经预滤波的信号进行多小波分解；
%
%  Input:  
%    maxlevel    integer, number of levels of decomposition
%    flt         string of characters, name of the filter;
%                for possible names and short descriptions see coef.m
%    fp          r by n real array, preprocessed input data;
%                r is the number of scaling functions,   
%                n must be of the form integer*2^maxlevel
%
%  Output:
%    transf      r by n real array, multiscaling and multiwavelet coefficients;
%                transf is organized as follows:   
%                columns 1 to n/2^maxlevel -- multiscaling coefficients
%                columns n/2^maxlevel +1 to n/2^(maxlevel-1) -- coarsest 
%                   multiwavelet coefficients 
%                            ........................  
%                columns n/2+1 to n -- finest multiwavelet coefficient
%
%  Example of Usage:
%   transf=dec1D_pe(fp,'cl',5) 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

n=length(fp(1,:));

[L,H]=coef(flt);
[nf,ns]=size(L);
[nf,nw]=size(H);
ns=ns/nf; 
nw=nw/nf;
if ns < nw
  big=nw;
else
  big=ns;
end
if fix(big/2) ~= big/2
  big=(big+1)/2;
else
  big=big/2;
end

i=1;
nss=1;
while L(:,nf*(i-1)+1:nf*i)==zeros(nf)
  nss=nss+1;
  i=i+1;
end
i=1;
nws=1;
while H(:,nf*(i-1)+1:nf*i)==zeros(nf)
  nws=nws+1;
  i=i+1;
end

transf=fp;
for lv=1:maxlevel
  fp=zeros(nf,n);
  fp=transf(:,1:n);
  tr=zeros(nf,n);
  for i=1:n/2-big
    for j=nss:ns
      tr(:,i)=tr(:,i)+L(:,nf*(j-1)+1:nf*j)*fp(:,2*i-2+j);
    end
    for j=nws:nw
      tr(:,i+n/2)=tr(:,i+n/2)+H(:,nf*(j-1)+1:nf*j)*fp(:,2*i-2+j);
    end
  end

  for i=n/2-big+1:n/2
    for j=nss:ns
      if 2*i-2+j>n
        tr(:,i)=tr(:,i)+L(:,nf*(j-1)+1:nf*j)*fp(:,2*i-2+j-n);
      else
        tr(:,i)=tr(:,i)+L(:,nf*(j-1)+1:nf*j)*fp(:,2*i-2+j);
      end
    end
    for j=nws:nw
      if 2*i-2+j>n
        tr(:,i+n/2)=tr(:,i+n/2)+H(:,nf*(j-1)+1:nf*j)*fp(:,2*i-2+j-n);
      else
        tr(:,i+n/2)=tr(:,i+n/2)+H(:,nf*(j-1)+1:nf*j)*fp(:,2*i-2+j);
      end
    end
  end
  transf(:,1:n)=tr;
  n=n/2;
end






