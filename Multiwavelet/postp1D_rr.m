function fhat=postp1D_rr(fhatp,flt)
%fhat=postp1D_rr(fhatp,flt)
%
% This function performs oversampled ("repeated row")  postprocessing of
% reconstructed preprocessed 1-dimensional signal. For details see [SW].
%
% 对重构后的信号进行后滤波；这个函数与prep1D_rr.m共同使用；
%
%  Input:                                    
%    flt        string of characters, name of the filter which was used; 
%               for possible names and short descriptions see coef.m
%    fhatp      r by n real array, reconstructed preprocessed signal;
%               each row corresponds to a scaling function  //这是重构后的信号；
%
%  Output: 
%    fhat       1 by n real array, reconstructed signal
%
%  Example of Usage:
%    fhat=postp1D_rr(fhatp,'ghm')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);
[nf,ns]=size(L);
ns=ns/nf;

H0=zeros(nf);
for i=1:ns
  H0=H0+L(:,nf*(i-1)+1:nf*i);
end
[V,D]=eig(H0'/(sqrt(2)));
fl=0;
for i=1:nf
  if abs(D(i,i)-1.)<0.000001
    fl=i;
  end
end
a0=V(:,fl);
big=0;
fl=0;
for i=1:nf
  if abs(a0(i))>abs(big)
    big=a0(i);
    fl=i;
  end
end

fhat=fhatp(fl,:);





