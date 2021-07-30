function transfth=thresh_vec2D(transf,maxlevel,thresh,covname,threshtype)
%transfth=thresh_vec2D(transf,maxlevel,thresh,covname,threshtype)
%
%  This function performs vector thresholding with decorrelation of 
%  2-dimensional multiwavelet transform. Scaling coefficients stay untouched. 
%  For description of the algorithm see [DS], [SW].
% 
%  Input:  
%    threshtype  string of characters, either 'hard' or 'soft';
%                  threshtype='hard' invokes hard thresholding,
%                  threshtype='soft' invokes soft thresholding;
%                for details see [DS]
%    covname     string of characters, name of covariance matrix used for 
%                decorrelation; for possible names and short descriptions 
%                see  coef_dcov2D.m
%    thresh      real, the value of the threshold
%    maxlevel    integer, the number of levels of wavelet decomposition 
%                in transf, maxlevel < 6
%    transf      n by n real array, wavelet transform of a signal;
%                n must be of the form integer*2^maxlevel;
%                for structure of transf see [SW]
%
%  Output: 
%    transfth    n by n real array, thresholded multiwavelet coefficients
%
%  Example of Usage:
%    transfth=thresh_vec2D(transf,3,0.75,'ghmap','hard')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

n0=length(transf(1,:));
[HH,HL,LH]=coef_dcov2D(covname);
r=length(HH(:,1));
nf=round(sqrt(r));

for i=1:maxlevel,
  hh=inv(HH(:,(i-1)*r+1:i*r));
  hl=inv(HL(:,(i-1)*r+1:i*r));
  lh=inv(LH(:,(i-1)*r+1:i*r));
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
      val=sqrt(tj'*hh*tj);
      if val<thresh
        tj=zeros(r,1);
        for m1=1:nf
          for m2=1:nf
            transf(j+(m1-1)*ste,k+(m2-1)*ste)=0;
          end
        end
      end
      if strcmp(threshtype,'soft'),
        tj=tj*(val-thresh)/val;
        for m1=1:nf
          for m2=1:nf
            transf(j+(m1-1)*ste,k+(m2-1)*ste)=tj((m1-1)*nf+m2);
          end
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
      val=sqrt(tj'*lh*tj);
      if val<thresh
        tj=zeros(r,1);
        for m1=1:nf
          for m2=1:nf
            transf(j+(m1-1)*ste,k+(m2-1)*ste)=0;
          end
        end
      end
      if strcmp(threshtype,'soft'),
        tj=tj*(val-thresh)/val;
        for m1=1:nf
          for m2=1:nf
            transf(j+(m1-1)*ste,k+(m2-1)*ste)=tj((m1-1)*nf+m2);
          end
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
      val=sqrt(tj'*hl*tj);
      if val<thresh
        tj=zeros(r,1);
        for m1=1:nf
          for m2=1:nf
            transf(j+(m1-1)*ste,k+(m2-1)*ste)=0;
          end
        end
      end
      if strcmp(threshtype,'soft'),
        tj=tj*(val-thresh)/val;
        for m1=1:nf
          for m2=1:nf
            transf(j+(m1-1)*ste,k+(m2-1)*ste)=tj((m1-1)*nf+m2);
          end
        end
      end      
    end
  end
end

transfth=transf;


