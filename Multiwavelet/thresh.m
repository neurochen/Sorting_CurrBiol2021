function transfth=thresh(transf,maxlevel,thresh,threshtype)
%transfth=thresh(transf,maxlevel,thresh,threshtype)
%
%  This function performs soft or hard thresholding of 1 or 2-dimensional 
%  wavelet transform. Scaling coefficients are not thresholded.
%  For description of soft and hard thresholding methods see [DJ].
%
%  Input:
%    threshtype  string of characters, either 'hard' or 'soft', 
%                  threshtype='hard' invokes hard thresholding,  
%                  threshtype='soft' invokes soft thresholding,    
%    thresh      real, value of the threshold 
%    maxlevel    integer, number of levels of wavelet decomposition in transf
%    transf      m by n reall array, wavelet transform of a signal;
%                m is either the number of scaling functions or m=n,
%                n must be of the form integer*2^maxlevel;
%                for structure of transf see [SW]
%
%  Output:
%    transfth    m by n real array, thresholded wavelet coefficients
%
% Example od Usage:
%   transfth=thresh(transf,3,0.75,'hard')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela


[n01,n02]=size(transf);
transfth=transf;
n2=n02/2^maxlevel;
n1=1;
if n01==n02,
  n1=n2;
end

for i=n1+1:n01,
  for j=n2+1:n02,
    if abs(transf(i,j))<thresh, 
      transfth(i,j)=0.;
    end
    if strcmp(threshtype,'soft'),
      if transf(i,j)<-thresh
        transfth(i,j)= transf(i,j)+thresh;
      end
      if  transf(i,j)>thresh,
        transfth(i,j)= transf(i,j)-thresh;
      end
    end
  end
end

if n01==n02,
  for i=1:n1,
    for j=n2+1:n02,
      if abs(transf(i,j))<thresh, 
        transfth(i,j)=0.;
      end
      if strcmp(threshtype,'soft'),
        if transf(i,j)<-thresh
          transfth(i,j)= transf(i,j)+thresh;
        end
        if  transf(i,j)>thresh,
          transfth(i,j)= transf(i,j)-thresh;
        end
      end
    end
  end
  for i=n1+1:n01,
    for j=1:n2,
      if abs(transf(i,j))<thresh, 
        transfth(i,j)=0.;
      end
      if strcmp(threshtype,'soft'),
        if transf(i,j)<-thresh
          transfth(i,j)= transf(i,j)+thresh;
        end
        if  transf(i,j)>thresh,
          transfth(i,j)= transf(i,j)-thresh;
        end
      end
    end
  end
end




