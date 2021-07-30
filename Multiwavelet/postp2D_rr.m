function fhat=postp2D_rr(fhatp,flt)
%fhat=postp2D_rr(fhatp,flt)
%
%  This function performs oversampled ("repeated row")  postprocessing of 
%  reconstructed preprocessed 2-dimensional signal. For details see [SW].
%
%  Input:   
%    flt        string of characters, name of the filter which was used; 
%               for possible names and short descriptions see coef.m 
%    fhatp      n by n real array, reconstructed preprocessed signal;
%               for structure of fhatp see [SW] 
%
%  Output:
%    fhat       n by n real array, reconstructed signal
%
%  Example of Usage:
%    fhat=postp2D_rr(fhatp,'ghm')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);
[n1,n2]=size(fhatp);
[nf,ns]=size(L);
ns=ns/nf;
n1=n1/nf;
n2=n2/nf;
fhat=zeros(n1,n2);

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

aa=zeros(nf,n1);
ff=zeros(n1,nf*n2);
for i=1:n2*nf
  for j=1:nf
    aa(j,:)=fhatp((j-1)*n1+1:j*n1,i)';
  end
  ff(:,i)=(aa(fl,:))';
end

aa=zeros(nf,n2);
for i=1:n1
  for j=1:nf
    aa(j,:)=ff(i,(j-1)*n2+1:j*n2);
  end
  fhat(i,:)=aa(fl,:);
end