function transf=dec2D_pe(fp,flt,maxlevel)
%transf=dec2D_pe(fp,flt,maxlevel)
%
%  This function computes discrete multiwavelet transform of a preprocessed 
%  2-dimensional signal. Boundaries are handled by periodic extension.
%  For description of the algorithm see [SW].
%
%  Input: 
%    maxlevel    integer, number of levels of decomposition
%    flt         string of characters, name of the filter;
%                for possible names and short descriptions see coef.m
%    fp          n by n real array, preprocessed input data;
%                n must be of the form integer*2^maxlevel;
%                for details of organization of fp see [SW].
%                                                   
%  Output:  
%    transf      n by n real array, multiscaling and multiwavelet coefficients;
%                for details of organization of transf see [SW] 
%
%  Example of Usage:
%   transf=dec2D_pe(fp,'bighm2',3)  

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[n1,n2]=size(fp);

[L,H]=coef(flt);
[nf,ns]=size(L);

transf=fp;
for lv=1:maxlevel
  nh=n2/nf;
  aa=zeros(nf,nh);
  for l=1:n1
    for j=1:nf
      aa(j,:)=transf(l,(j-1)*nh+1:j*nh);
    end
    aa=dec1D_pe(aa,flt,1);
    for j=1:nf
      transf(l,(j-1)*nh/2+1:j*nh/2)=aa(j,1:nh/2);
      transf(l,n2/2+(j-1)*nh/2+1:n2/2+j*nh/2)=aa(j,nh/2+1:nh);    
    end
  end
  nh=n1/nf;
  aa=zeros(nf,nh);
  for l=1:n2
    for j=1:nf
      aa(j,:)=transf((j-1)*nh+1:j*nh,l)';
    end
    aa=dec1D_pe(aa,flt,1);
    for j=1:nf
      transf((j-1)*nh/2+1:j*nh/2,l)=aa(j,1:nh/2)';
      transf(n1/2+(j-1)*nh/2+1:n1/2+j*nh/2,l)=aa(j,nh/2+1:nh)';    
    end
  end
  n1=n1/2;
  n2=n2/2;
end