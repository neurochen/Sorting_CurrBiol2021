function transfth=thresh_decor2D(transf,maxlevel,thresh,covname,threshtype)
%transfth=thresh_decor2D(transf,maxlevel,thresh,covname,threshtype)
%
%  This function performs scalar thresholding with decorrelation of 
%  2-dimensional wavelet transform. Scaling coefficients stay untouched.
%  Wavelet coefficients are decorrelated using diagonal blocks of
%  the covariance matrix of the transform. After that all coefficients are 
%  thresholded and correlated again.  
%
%  Input: 
%    threshtype  string of characters, either 'hard' or 'soft';
%                  threshtype='hard' invokes hard thresholding, 
%                  threshtype='soft' invokes soft thresholding,
%                for details see [DJ] in VSMWP.README
%    covname     string of characters, name of covariance matrix used for 
%                decorrelation; for possible names and short descriptions
%                see coef_dcov1D.m
%    thresh      real, value of the threshold 
%    maxlevel    integer, number of levels of wavelet decomposition in transf;
%                maxleve < 6 
%    transf      n by n real array, wavelet transform of a signal;
%                n must be of the form integer*2^maxlevel;
%                for structure of transf see [SW]
%
%  Output: 
%    transfth    n by n real array, thresholded wavelet coefficients
%
%  Example of Usage:
%    transfth=thresh_decor2D(transf,5,0.75,'ghmap','hard')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

n0=length(transf(1,:));
[HH,HL,LH]=coef_dcov2D(covname);
r=length(HH(:,1));
nf=round(sqrt(r));

for i=1:maxlevel,
  hh=inv(sqrtm(HH(:,(i-1)*r+1:i*r)));
  ihh=sqrtm(HH(:,(i-1)*r+1:i*r));
  hl=inv(sqrtm(HL(:,(i-1)*r+1:i*r)));
  ihl=sqrtm(HL(:,(i-1)*r+1:i*r));
  lh=inv(sqrtm(LH(:,(i-1)*r+1:i*r)));
  ilh=sqrtm(LH(:,(i-1)*r+1:i*r));
  nb=n0-(2^i-1)*n0/2^i;
  ste=round(n0/(nf*2^i));

  for j=nb+1:nb+ste,
    for k=nb+1:nb+ste,
      tj=zeros(r,1);
      for m1=1:nf
        for m2=1:nf
          tj((m1-1)*nf+m2)=transf(j+(m1-1)*ste,k+(m2-1)*ste);
        end
      end
      tj=hh*tj;
      for kk=1:r
        if abs(tj(kk))<thresh, 
          tj(kk)=0.;
        end
        if strcmp(threshtype,'soft'),
          if tj(kk)<-thresh
            tj(kk)=tj(kk)+thresh;
          end
          if tj(kk)>thresh,
            tj(kk)=tj(kk)-thresh;
          end
        end
      end
      tj=ihh*tj;
      for m1=1:nf
        for m2=1:nf
          transf(j+(m1-1)*ste,k+(m2-1)*ste)=tj((m1-1)*nf+m2);
        end
      end
    end
  end

  for j=1:ste
    for k=nb+1:nb+ste,
      tj=zeros(r,1);
      for m1=1:nf
        for m2=1:nf
          tj((m1-1)*nf+m2)=transf(j+(m1-1)*ste,k+(m2-1)*ste);
        end
      end
      tj=lh*tj;
      for kk=1:r
        if abs(tj(kk))<thresh, 
          tj(kk)=0.;
        end
        if strcmp(threshtype,'soft'),
          if tj(kk)<-thresh
            tj(kk)=tj(kk)+thresh;
          end
          if tj(kk)>thresh,
            tj(kk)=tj(kk)-thresh;
          end
        end
      end
      tj=ilh*tj;
      for m1=1:nf
        for m2=1:nf
          transf(j+(m1-1)*ste,k+(m2-1)*ste)=tj((m1-1)*nf+m2);
        end
      end
    end
  end

  for j=nb+1:nb+ste,
    for k=1:ste,
      tj=zeros(r,1);
      for m1=1:nf
        for m2=1:nf
          tj((m1-1)*nf+m2)=transf(j+(m1-1)*ste,k+(m2-1)*ste);
        end
      end
      tj=hl*tj;
      for kk=1:r
        if abs(tj(kk))<thresh, 
          tj(kk)=0.;
        end
        if strcmp(threshtype,'soft'),
          if tj(kk)<-thresh
            tj(kk)=tj(kk)+thresh;
          end
          if tj(kk)>thresh,
            tj(kk)=tj(kk)-thresh;
          end
        end
      end
      tj=ihl*tj;
      for m1=1:nf
        for m2=1:nf
          transf(j+(m1-1)*ste,k+(m2-1)*ste)=tj((m1-1)*nf+m2);
        end
      end
    end
  end
end

transfth=transf;

