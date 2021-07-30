function [rmsecomp,fhat]=simpcomp2D(fname,f,dflt,rflt,prepr,nlarge)
%[rmsecomp,fhat]=simpcomp2D(fname,f,dflt,rflt,prepr,nlarge)
%
%  This function compresses 2-dimensional signal by retaining given number of 
%  largest coefficients in its wavelet transform. Compressed transform is
%  reconstructed and the result is compared to the initial signal. 
%  For details see  [SW1].
%
%  [SW1] V. Strela and A. T. Walden, "Orthogonal and biorthogonal multiwavelets 
%        for signal denoising and image compression", SPIE Proc. 3391 
%        AeroSense 98, Orlando, Florida, April 1998.
%
%
%  Input:                                                   
%    nlarge      integer, number of wavelet coefficients to retain
%    prepr       string of characters, type of preprocessing; 
%                allowed are:
%                  any prefilter names from coef_prep.m, they imply critically 
%                  sampled preprocessing and should be used together with 
%                  appropriate multifilter,
%                  prepr='rr'     corresponds to oversampled ("repeated row")
%                                 preprocessing and  can be used with any 
%                                 multifilter,
%                  prepr='scalar' indicates that the filter is scalar and no
%                                 preprocessing is needed, 
%                for details see [SW] in VSMWP.README
%    rflt        string of characters, name of the reconstruction filter bank;
%                for admissible names see coef.m  //重建滤波器
%    dflt        string of characters, name of the decomposition filter bank;
%                for admissible names see coef.m  //分解滤波器
%    f           n by m real array, input signal, n and m must be of 
%                the form 2^k, k > 4 is integer   //输入信号
%    name        string of characters, name of the input data (used only
%                to create the name of the output PostScript file)  //输入信号的名称      
%
%  Output: 
%    fhat        n by m real array, reconstructed signal  //重建信号
%    rmsecomp    real, root mean square error of the reconstructed signal
%
%  Initial and reconstructed signals are plotted, reconstructed signal is 
%  saved in a PostScript file.
%
%  Example of Usage:
%  [r,A]=simpcomp2D('my_image',f,'bih52s','bih32s','bih5ap',1024);

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

[n1,n2]=size(f);
n0=n1;
if n2<n1
  n0=n2;
end
maxlevel=round(log(n0)/log(2))-4;

signal=fname
total_number_of_coefficients=n1*n2
number_of_coefficients_to_keep=nlarge
transform=strcat(dflt,'_',rflt,'_',prepr)

figure(1);
colormap(gray(256))
image(f)
axis('square')
axis('off')
title('Original Image')

if strcmp(prepr,'rr')
  fp=prep2D_rr(f,dflt);
else
  if strcmp(prepr,'scalar')
    fp=f;
  else 
    fp=prep2D_appe(f,prepr);
    maxlevel=maxlevel;
  end
end

transf=dec2D_pe(fp,dflt,maxlevel);

transfcomp=keepnlarge(transf,nlarge);

fhatp=rec2D_pe(transfcomp,rflt,maxlevel);

if strcmp(prepr,'rr')
  fhat=postp2D_rr(fhatp,dflt);
else
  if strcmp(prepr,'scalar')
    fhat=fhatp;
  else 
    fhat=postp2D_appe(fhatp,prepr);
  end
end

rmsecomp=norm(f-fhat,'fro')/sqrt(n1*n2)
figure(2);
colormap(gray(256))
image(fhat)
axis('square')
axis('off')
title(strcat('Compressed Image: ',num2str(nlarge),' largest coefficients kept.'))

%psname=strcat('comp_',fname,'_',dflt,prepr,'_',num2str(nlarge),'.ps');
%print(psname);

%INFO=strcat('Compressed image is saved in --  ','   ',psname)

