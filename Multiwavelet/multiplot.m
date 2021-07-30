function mes=multiplot(flt)
%mes=multiplot(flt)
%
%  This routine plots scaling and wavelet functions corresponding to the given
%  multifilter bank. 
%
%  Input:
%    flt         string of characters, name of the filter bank;
%                for possible names and short descriptions see coef.m
%
%  Output:       2*n plots: n scaling functions and n wavelets
%
%  Example of Usage:
%    multiplot('ghm')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

lmax=10;
n=2^lmax;
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
mes=strcat('wait until all_',num2str(2*nf),' functions are plotted...');

lmin=ceil(log(big)/log(2));
maxlevel=lmax-lmin;

dum=ones(1,nf);
for i=1:nf
  aa=zeros(nf,n);
  aa(i,2^(lmin-1)-1)=1;
  f=rec1D_pe(aa,flt,maxlevel);
  figure(2*i-1)
  plot((2^lmin)*[0:n-1]/n,dum*f*2^(maxlevel/2));
  axis('tight')
%  title(strcat('scaling function #',num2str(i)))
  aa=zeros(nf,n);
  aa(i,3*2^(lmin-1)-1)=1;
  f=rec1D_pe(aa,flt,maxlevel);
  figure(2*i)
  plot((2^lmin)*[0:n-1]/n,dum*f*2^(maxlevel/2));
  axis('tight')
 % title(strcat('wavelet function #',num2str(i))) 
end 






