function fp=prep2D_rr(f,flt)
%fp=prep2D_rr(f,flt) 
%
%  This function performs oversampled ("repeated row") preprocessing of 
%  given 2-dimensional signal. For description of the algorithm see [SW].
%
%  Input:
%    flt        string of characters, name of the filter to be used;
%               for possible names and short descriptions see coef.m 
%    f          n by n real array, input signal;//2-D输入信号
%
%    n          n must be of the form integer*2^maxlevel,//信号的长度
%    maxlevel   maxlevel is the number of levels of multiwavelet decomposition //多小波分解的层数
%
%  Output: 
%    fp         n by n real array, preprocessed signal;
%               for structure of fp see [SW]
%
%  Example of Usage
%    fp=prep2D_rr(f,'cl') 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);
[n1,n2]=size(f);
[nf,ns]=size(L);
ns=ns/nf;
fp=zeros(nf*n1,nf*n2);

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

for i=1:n1
  for j=1:nf
    fp(i,(j-1)*n2+1:j*n2)=a0(j)*f(i,:);
  end
end

aa=fp(1:n1,:);
for i=1:n2*nf
  for j=1:nf
    fp((j-1)*n1+1:j*n1,i)=a0(j)*aa(:,i);
  end
end




