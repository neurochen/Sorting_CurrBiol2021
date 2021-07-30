function fp=prep1D_rr(f,flt)
%fp=prep1D_rr(f,flt) 
%
%  This function performs oversampled ("repeated row") preprocessing of given 
%  1-dimensional signal. For description of the algorithm see [SW]
%
%  rr 这是重复行预滤波；这个函数与post1D_rr.m共同使用；
%
%  Input:
%    flt        string of characters, name of the filter to be used; 
%               for possible names and short descriptions see coef.m 
%    f          1 by n real array, input signal;//1-D输入信号
%
%    n          n must be of the form integer*2^maxlevel,//输入信号的长度
%    maxlevel   maxlevel is the number of levels of multiwavelet decomposition //多小波分解的层数
%
%  Output: 
%    fp         r by n real array, preprocessed signal; 
%               each row corresponds to a scaling function 
%
%  Example of Usage:
%    fp=prep1D_rr(f,'ghm') 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);

n=length(f);
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
for i=1:nf
  if abs(a0(i))>abs(big)
    big=a0(i);
  end
end
for i=1:nf
  a0(i)=a0(i)/big;
end

fp=zeros(nf,n);
for i=1:nf
  fp(i,:)=a0(i)*f;
end