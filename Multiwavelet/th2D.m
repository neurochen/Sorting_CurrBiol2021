function [rmsedenoise,fhat]=th2D(fname,f,dflt,rflt,prepr,threshmeth,threshtype,SNR)
%[rmsedenoise,fhat]=th2D(fname,f,dflt,rflt,prepr,threshmeth,threshtype,SNR)
%
%  This function adds Gaussian white noise to 2-dimensional signal,  
%  performs denoising by given thresholding method and compares the result 
%  to the initial signal. For details see [SW].
%
%  Input:                                                   
%    SNR         real, signal to noise ratio (ratio of the standard deviation
%                of the signal to the standard deviation of the noise)
%    threshtype  string of characters, either 'hard' or 'soft', 
%                  threshtype='hard' invokes hard thresholding,   
%                  threshtype='soft' invokes soft thresholding,   
%                for details see [DJ] 
%    threshmeth  string of characters, either 'scalar' or 'vector' or 'decor',
%                  threshmeth='scalar' adjusts Donoho threshold (see [DJ])
%                                      by mean variance of the transform, 
%                                      and starts scalar thresholding, see [SW]
%                  threshmeth='vector' starts vector thresholding with 
%                                      decorrelation, see [DS], [SW]
%
%
%    [DS]  T. R. Downie and B. W. Silverman, "The discrete multiple wavelet
%      transform and thresholding methods", IEEE Trans. on SP, to appear
%                  threshmeth='decor'  starts scalar thresholding with decorrelation 
%
%    prepr       string of characters, type of preprocessing; 
%                allowed are:
%                  prefilter names from coef_prep.m, they correspond
%                  to critically sampled preprocessing and should be used 
%                  together with appropriate multifilter,
%                  prepr='rr'      corresponds to "repeated row" (oversampled)
%                                  preprocessing and  can be used with any 
%                                  multifilter,
%                  prepr= 'scalar' indicates that the filter is scalar and no
%                                  preprocessing is needed, 
%                for details see [SW]
%    rflt        string of characters, name of the reconstruction filter bank;
%                for admissible names see coef.m
%    dflt        string of characters, name of the decomposition filter bank;
%                for admissible names see coef.m
%    f           n by n real array, input signal; 
%                n must be of the form 2^k, 4 < k < 11 is integer 
%    name        string of characters, name of the input data (used only
%                to create a name of the output PostScript file)        
%
%  Output: 
%    fhat        n by n real array, denoised signal
%    rmsedenoise real, root mean square error of the denoised signal
%
%  Input and denoised signals are plotted, denoised signal is saved in a 
%  PostScript file.
%  WARNING: Not all combinations of filters, prefilters and types of 
%  thresholding can be used!
%
%  Example of Usage:
%    [r,A]=th2D('my_image',f,'bih52s','bih32s','bih5ap','decor','hard',2)

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

signal=fname
transform=strcat(dflt,'_',rflt,'_',prepr)
thresholding=strcat(threshmeth,'_',threshtype)
SNR

figure(1);
colormap(gray(256))
image(f)
axis('square')
axis('off')
%title('Original Image')

n0=length(f(1,:));
maxlevel=round(log(n0)/log(2))-4;

sig=norm(f-mean(mean(f)),'fro')/n0;
sig=sig/SNR;

randn('state',0);
noise=randn(n0);
fn=f+sig*noise;
noise_rmse=norm(f-fn,'fro')/n0

figure(2);
colormap(gray(256))
image(fn)
axis('square')
axis('off')
%title(strcat('Noisy Image, SNR=',num2str(SNR)))
%psname=strcat('th_',fname,'_noisy_',num2str(SNR),'.ps');
%print(psname);
%INFO=strcat('Noisy image is saved in --  ','   ',psname)

thr=sqrt(2*log(n0))*sig;
if strcmp(prepr,'rr')
  thr=sqrt(2*log(2*n0))*sig;
end

pflt=strcat(dflt,prepr);

if strcmp(threshmeth,'scalar') 
  thadj=coef_thadj(pflt);
  thr=thadj(round(log(n0)/log(2))-4,maxlevel)^2*thr;
end

if strcmp(prepr,'rr')
  fp=prep2D_rr(fn,dflt);
else
  if strcmp(prepr,'scalar')
    fp=fn;
    pflt=dflt;
  else 
    fp=prep2D_appe(fn,prepr);
    maxlevel=maxlevel-1;
    pflt=prepr;
  end
end

transf=dec2D_pe(fp,dflt,maxlevel);

if strcmp(threshmeth,'scalar')
  transftr=thresh(transf,maxlevel,thr,threshtype);
elseif strcmp(threshmeth,'vector')
  transftr=thresh_vec2D(transf,maxlevel,thr,pflt,threshtype);
elseif strcmp(threshmeth,'decor')
  transftr=thresh_decor2D(transf,maxlevel,thr,pflt,threshtype);
end

fhatp=rec2D_pe(transftr,rflt,maxlevel);

if strcmp(prepr,'rr')
  fhat=postp2D_rr(fhatp,dflt);
else
  if strcmp(prepr,'scalar')
    fhat=fhatp;
  else 
    fhat=postp2D_appe(fhatp,prepr);
  end
end

rmsedenoise=norm(f-fhat,'fro')/n0
figure(3);
colormap(gray(256))
image(fhat)
axis('square')
axis('off')
%title(strcat('Denoised Image, rmse=',num2str(rmsedenoise)))
%psname=strcat('th_',fname,'_',dflt,prepr,'_',threshmeth,'_',threshtype,'_',num2str(SNR),'.ps');
%print(psname);

%INFO=strcat('Denoised image is saved in --  ','   ',psname)






