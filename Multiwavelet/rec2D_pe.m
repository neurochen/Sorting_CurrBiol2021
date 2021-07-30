function fhatp=rec2D_pe(transf,flt,maxlevel)
%fhatp=rec2D_pe(transf,flt,maxlevel)
%
%  This function reconstructs preprocessed 2-dimensional signal from its 
%  discrete multiwavelet transform. Boundaries are handled by periodic 
%  extension. For description of the algorithm see [SW]
%
%  Input:
%    maxlevel    integer, number of levels of decomposition 
%    flt         string of characters, name of the filter bank;
%                for possible names and short descriptions see coef.m
%    transf      n by n real array, multiscaling and multiwavelet coefficients;
%                for structure of transf see [SW] 
%                n must be of the form integer*2^maxlevel
%
%  Output: 
%    fhatp       n by n real array, reconstructed preprocessed signal;
%                for structure of fhatp see [SW]
%
%  Example of Usage:
%    fhatp=rec2D_pe(transf,'cl',3)

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[n1,n2]=size(transf);
n1=n1/2^(maxlevel-1);
n2=n2/2^(maxlevel-1);

[L,H]=coef(flt);
[nf,ns]=size(L);

fhatp=transf;
for lv=1:maxlevel
  nh=n1/nf;
  aa=zeros(nf,nh);
  for l=1:n2
    for j=1:nf
      aa(j,1:nh/2)=fhatp((j-1)*nh/2+1:j*nh/2,l)';
      aa(j,nh/2+1:nh)=fhatp(n1/2+(j-1)*nh/2+1:n1/2+j*nh/2,l)';    
    end
    aa=rec1D_pe(aa,flt,1);
    for j=1:nf
      fhatp((j-1)*nh+1:j*nh,l)=aa(j,:)';
    end
  end
  nh=n2/nf;
  aa=zeros(nf,nh);
  for l=1:n1
    for j=1:nf
      aa(j,1:nh/2)=fhatp(l,(j-1)*nh/2+1:j*nh/2);
      aa(j,nh/2+1:nh)=fhatp(l,n2/2+(j-1)*nh/2+1:n2/2+j*nh/2);    
    end
    aa=rec1D_pe(aa,flt,1);
    for j=1:nf
      fhatp(l,(j-1)*nh+1:j*nh)=aa(j,:);
    end
  end
  n1=n1*2;
  n2=n2*2;
end
