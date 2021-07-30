function T=transop(flt)
%T=transop(flt) 
%
%  This function computes the matrix of the transition operator corresponding 
%  to a given multiscaling function.
%  If c is the largest eigenvalue of T which is either non-dyadic or multiple 
%  dyadic, then Sobolev regularity of the multiscaling function is higher 
%  than -log(c)/log(4).      
%  For theory and better estimates of Sobolev regularity see [J].
% 
%  Input:
%    flt         string of characters, name of the filter;
%                for possible names and short descriptions see coef.m
%
%  Output: 
%    T           r^2*(2*l-1) by r^2*(2*l-1) real array, matrix of the 
%                transition operator;
%                r is the number of scaling functions,
%                l is the number of scaling coefficients  
%
%  Example of Usage:
%    T=transop('ghm')

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[L,H]=coef(flt);
r=length(L(:,1));
n=length(L(1,:))/r;
L2=zeros(r*r,(2*n-1)*r*r);

for i=1:n,
   LL=kron(L(:,(i-1)*r+1:i*r),L(:,(i-1)*r+1:i*r));
   L2(:,(n-1)*r*r+1:n*r*r)=L2(:,(n-1)*r*r+1:n*r*r)+LL;
end
for i=1:n-1,
  for j=1:i,
   LL=kron(L(:,(n-i+j-1)*r+1:(n-i+j)*r),L(:,(j-1)*r+1:j*r));
   L2(:,(i-1)*r*r+1:i*r*r)=L2(:,(i-1)*r*r+1:i*r*r)+LL;
   LL=kron(L(:,(j-1)*r+1:j*r),L(:,(n-i+j-1)*r+1:(n-i+j)*r));
   L2(:,(2*n-i-1)*r*r+1:(2*n-i)*r*r)=L2(:,(2*n-i-1)*r*r+1:(2*n-i)*r*r)+LL;
  end,
end

r1=length(L2(:,1));
n1=length(L2(1,:))/r1;

T1=zeros(n1*r1,n1*r1); 
for i=1:n1, 
  for j=1:n1,
    if (0<2*i-j)&(2*i-j<=n1) 
      T1((i-1)*r1+1:i*r1,(j-1)*r1+1:j*r1)=L2(1:r1,(2*i-j-1)*r1+1:(2*i-j)*r1);
    end
  end,
end

T=T1(r*r+1:(2*n-2)*r*r,r*r+1:(2*n-2)*r*r);

