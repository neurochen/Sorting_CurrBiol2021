function fhatp=rec1D_pe(transf,flt,maxlevel)
%fhatp=rec1D_pe(transf,flt,maxlevel)
%
%  This function reconstructs preprocessed 1-dimensional signal from its 
%  discrete multiwavelet transform. Boundaries are handled by periodic 
%  extension. For description of the algorithm see [SW].
%
%  对已经进行多小波分解的信号进行重构
%
%  Input:  
%    maxlevel    integer, number of levels of decomposition 
%    flt         string of characters, name of the filter bank;
%                for possible names and short descriptions see coef.m
%    transf      r by n real array, multiscaling and multiwavelet coefficients;
%                transf is organized as follows:
%                columns 1 to n/2^maxlevel -- multiscaling coefficients
%                columns n/2^maxlevel +1 to n/2^(maxlevel-1) -- coarsest 
%                        multiwavelet coefficients
%                            ........................    
%                columns n/2+1 to n -- finest multiwavelet coefficients;
%                r is the number of scaling functions,  
%                n must be of the form integer*2^maxlevel       
%
%  Output:  
%    fhatp       r by n real array, reconstructed preprocessed signal
%
%  Example of Usage:
%    fhatp=rec1D_pe(transf,'cl',4)

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

n=length(transf(1,:))/2^(maxlevel-1);

[L,H]=coef(flt);
[nf,ns]=size(L);
[nf,nw]=size(H);
ns=ns/nf; 
nw=nw/nf;
if fix(ns/2) ~= ns/2
  nse=(ns+1)/2;
else
  nse=ns/2;
end
if fix(nw/2) ~= nw/2
  nwe=(nw+1)/2;
else
  nwe=nw/2;
end
nso=ns-nse;
nwo=nw-nwe;
if ns < nw
  big=nwe;
else
  big=nse;
end

i=1;
nses=1;
while L(:,nf*(2*i-2)+1:nf*(2*i-1))==zeros(nf)
  nses=nses+1;
  i=i+1;
end
i=1;
nsos=1;
while L(:,nf*(2*i-1)+1:nf*2*i)==zeros(nf)
  nsos=nsos+1;
  i=i+1;
end
i=1;
nwes=1;
while H(:,nf*(2*i-2)+1:nf*(2*i-1))==zeros(nf)
  nwes=nwes+1;
  i=i+1;
end
i=1;
nwos=1;
while H(:,nf*(2*i-1)+1:nf*2*i)==zeros(nf)
  nwos=nwos+1;
  i=i+1;
end


for lv=1:maxlevel
  fhatp=zeros(nf,n);
  for i=1:big-1
    for j=nses:nse
      if i-j < 0
        fhatp(:,2*i-1)=fhatp(:,2*i-1)+L(:,nf*(2*j-2)+1:nf*(2*j-1))'*transf(:,n/2+i-j+1);
      else      
        fhatp(:,2*i-1)=fhatp(:,2*i-1)+L(:,nf*(2*j-2)+1:nf*(2*j-1))'*transf(:,i-j+1);
      end
    end
    for j=nwes:nwe
      if i-j < 0
        fhatp(:,2*i-1)=fhatp(:,2*i-1)+H(:,nf*(2*j-2)+1:nf*(2*j-1))'*transf(:,n+i-j+1);
      else
        fhatp(:,2*i-1)=fhatp(:,2*i-1)+H(:,nf*(2*j-2)+1:nf*(2*j-1))'*transf(:,n/2+i-j+1);
      end
    end
    for j=nsos:nso
      if i-j < 0
        fhatp(:,2*i)=fhatp(:,2*i)+L(:,nf*(2*j-1)+1:nf*(2*j))'*transf(:,n/2+i-j+1);
      else
        fhatp(:,2*i)=fhatp(:,2*i)+L(:,nf*(2*j-1)+1:nf*(2*j))'*transf(:,i-j+1);
      end
    end
    for j=nwos:nwo
      if i-j < 0
        fhatp(:,2*i)=fhatp(:,2*i)+H(:,nf*(2*j-1)+1:nf*(2*j))'*transf(:,n+i-j+1);
      else
        fhatp(:,2*i)=fhatp(:,2*i)+H(:,nf*(2*j-1)+1:nf*(2*j))'*transf(:,n/2+i-j+1);
      end
    end
  end

  for i=big:n/2
    for j=nses:nse
      fhatp(:,2*i-1)=fhatp(:,2*i-1)+L(:,nf*(2*j-2)+1:nf*(2*j-1))'*transf(:,i-j+1);
    end
    for j=nwes:nwe
      fhatp(:,2*i-1)=fhatp(:,2*i-1)+H(:,nf*(2*j-2)+1:nf*(2*j-1))'*transf(:,n/2+i-j+1);
    end
    for j=nsos:nso
      fhatp(:,2*i)=fhatp(:,2*i)+L(:,nf*(2*j-1)+1:nf*(2*j))'*transf(:,i-j+1);
    end
    for j=nwos:nwo
      fhatp(:,2*i)=fhatp(:,2*i)+H(:,nf*(2*j-1)+1:nf*(2*j))'*transf(:,n/2+i-j+1);
    end
  end
  transf(:,1:n)=fhatp;
  n=2*n;
end




