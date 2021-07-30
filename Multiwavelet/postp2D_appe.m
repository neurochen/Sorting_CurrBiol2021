function fhat=postp2D_appe(fhatp,pflt)
%fhat=postp2D_appe(fhatp,pflt)
%
%  This function performs critically sampled postprocessing of reconstructed 
%  preprocessed 2-dimensional signal. Boundaries are handled by periodic 
%  extension. For details see [SW].
%
%  Input: 
%    pflt       string of characters, name of the prefilter;
%               for possible names and short descriptions see coef_prep.m
%    fhatp      n by n real array, reconstructed preprocessed signal; 
%               for structure of fhatp see [SW]
%
%  Output: 
%    fhat       n by n real array, reconstructed signal 
%
%  Example of Usage:
%    fhat=postp2D_appe(fhatp,'clap')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[PR,PO]=coef_prep(pflt);
[n1,n2]=size(fhatp);
fhat=zeros(n1,n2);
[nf,np]=size(PO);

prl=round(n1/nf);
aa=zeros(nf,prl);
for i=1:n2
  for j=1:nf
    aa(j,:)=fhatp((j-1)*prl+1:j*prl,i)';
  end
  fhat(:,i)=(postp1D_appe(aa,pflt))';
end

prl=round(n2/nf);
aa=zeros(nf,prl);
for i=1:n1
  for j=1:nf
    aa(j,:)=fhat(i,(j-1)*prl+1:j*prl);
  end
  fhat(i,:)=postp1D_appe(aa,pflt);
end