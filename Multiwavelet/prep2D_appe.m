function fp=prep2D_appe(f,pflt)
%fp=prep2D_appe(f,pflt)
%
%  This function performs critically sampled preprocessing of given 
%  2-dimensional signal. Boundaries are handled by periodic extension.
%  For description of the algorithm see [SW].
%
%  Input:    
%    pflt       string of characters, name of the prefilter;
%               for possible names and short descriptions see coef_prep.m
%    f          n by n real array, input signal;//2-D输入信号
%    n          n must be of the form r*integer*2^maxlevel,
%
%    r          r is the number of scaling functions, //尺度函数的个数
%    maxlevel   maxlevel is the number of levels of multiwavelet decomposition //多小波分解的层数
%
%  Output: 
%    fp         n by n real array, preprocessed data;
%               for structure of fp see [SW]
%
%  Example of Usage:
%    fp=prep2D_appe(f,'ghmap')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[PR,PO]=coef_prep(pflt);
[n1,n2]=size(f);
fp=zeros(n1,n2);
[nf,np]=size(PR);

prl=round(n2/nf);
aa=zeros(nf,prl);
for i=1:n1
  aa=prep1D_appe(f(i,:),pflt);
  for j=1:nf
    fp(i,(j-1)*prl+1:j*prl)=aa(j,:);
  end
end

prl=round(n1/nf);
aa=zeros(nf,prl);
for i=1:n2
  aa=prep1D_appe(fp(:,i)',pflt);
  for j=1:nf
    fp((j-1)*prl+1:j*prl,i)=aa(j,:)';
  end
end