function fp=prep1D_appe(f,pflt)
%fp=prep1D_appe(f,pflt)
%
%  This function performs critically sampled preprocessing of given 
%  1-dimensional signal. Boundaries are handled by periodic extension.
%  For description of the algorithm see [SW] //对于边界是周期扩张
%
%  appe: 这是逼近预滤波；这个函数与postp1D_appe.m共同使用；
%
%  Input:
%    pflt       string of characters, name of the prefilter;
%               for possible names and short descriptions see coef_prep.m
%    f          1 by n real array, input data; //1-D输入信号
%
%    n          n must be of the form r*integer*2^maxlevel, //输入信号的长度
%    r          r is the number of scaling functions,//尺度函数的个数
%    maxlevel   maxlevel is the number of levels of multiwavelet transform //多小波分解的层数
%
%  Output:  
%    fp         r by n/r real array, preprocessed signal;
%               each row corresponds to a scaling function 
%
%  Example of Usage:
%    fp=prep1D_appe(f,'clap')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[PR,PO]=coef_prep(pflt); %% 得到预滤波矩阵

n=length(f);

[nf,np]=size(PR);
np=np/nf;

fs=zeros(nf,round(n/nf));
for i=1:round(n/nf)
  for j=1:nf
    fs(j,i)=f(nf*(i-1)+j);
  end
end

fp=zeros(nf,round(n/nf));
for i=1:round(n/nf)-np+1
  for j=1:np
    fp(:,i)=fp(:,i)+PR(:,nf*(j-1)+1:nf*j)*fs(:,i+j-1);
  end
end

for i=round(n/nf)-np+2:round(n/nf)
  for j=1:np
    if i+j-1>round(n/nf)
      fp(:,i)=fp(:,i)+PR(:,nf*(j-1)+1:nf*j)*fs(:,i+j-1-round(n/nf));
    else
      fp(:,i)=fp(:,i)+PR(:,nf*(j-1)+1:nf*j)*fs(:,i+j-1);
    end
  end
end

