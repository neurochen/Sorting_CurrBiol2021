function transfcomp=keepnlarge(transf,nbigc)
%transfcomp=keepnlarge(transf,nbigc)
%
%  This function keeps given number of largest elements in an input array and
%  sets the rest of the array to zero. It can be used as the simplest
%  compression tool and for comparison of performance of different transforms. 
%  For comparison of several scalar and multiwavelet transforms using this 
%  function see [SW1]. 
%
%  Input:                                                 
%    nbigc       integer, number of largest coefficients to keep
%    transf      m by n real array, transform of a signal 
%
%  Output:   
%    transfcomp  m by n real array, nbigc largest coefficients are the same
%                as in transf and the rest is zero
%
%  Example of Usage:
%   transfcomp=keepnlarge(transf,1024) 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela


[n1,n2]=size(transf);
N=n1*n2;
transfcomp=zeros(n1,n2);
lon=zeros(1,N);

if n1>n2
  for i=1:n2
    lon((i-1)*n1+1:i*n1)=abs(transf(:,i));
  end
  [slon,ind]=sort(lon);
  for i=1:nbigc
    a=rem(ind(N-i+1)-1,n1);
    b=ceil(ind(N-i+1)/n1);
    transfcomp(a+1,b)=transf(a+1,b);
  end
else
  for i=1:n1
    lon((i-1)*n2+1:i*n2)=abs(transf(i,:));
  end
  [slon,ind]=sort(lon);
  for i=1:nbigc
    a=rem(ind(N-i+1)-1,n2);
    b=ceil(ind(N-i+1)/n2);
    transfcomp(b,a+1)=transf(b,a+1);
  end  
end
