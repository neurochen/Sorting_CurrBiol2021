% MWMP contains the following functions:
%
%   coef           coefficients of various multifilter banks(多重滤波器组系数)
%   coef_prep      coefficients of various multiwavelet pre/post filters(预滤波与后滤波)
%   coef_dcov1D    diagonal blocks of covariance(协方差) matrices of various 
%                  1-dimensional multiwavelet transforms
%   coef_dcov2D    diagonal blocks of covariance matrices of various 
%                  2-dimensional multiwavelet transforms
%   coef_thadj     mean variances of various multiwavelet transforms
%    
%   prep1D_appe    1-dimensional critically sampled multiwavelet preprocessing
%   prep1D_rr      1-dimensional oversampled multiwavelet preprocessing
%   dec1D_pe       1-dimensional multiwavelet decomposition(一维多重小波正变换)
%   rec1D_pe       1-dimensional multiwavelet reconstruction(一维多重小波逆变换)
%   postp1D_rr     1-dimensional oversampled multiwavelet postprocessing
%   postp1D_appe   1-dimensional critically sampled multiwavelet postprocessing
%   prep2D_appe    2-dimensional critically sampled multiwavelet preprocessing
%   prep2D_rr      2-dimensional oversampled multiwavelet preprocessing
%   dec2D_pe       2-dimensional multiwavelet decomposition
%   rec2D_pe       2-dimensional multiwavelet reconstruction
%   postp2D_rr     2-dimensional oversampled multiwavelet postprocessing
%   postp2D_appe   2-dimensional critically sampled multiwavelet postprocessing
%
%   keepnlarge     keeps given number of largest coefficients of a transform 
%   thresh         scalar thresholding of 1 or 2 dimensional multiwavelet 
%                  transform
%   thresh_decor1D scalar thresholding with decorrelation of 1-dimensional
%                  multiwavelet transform
%   thresh_vec1D   vector thresholding with decorrelation of 1-dimensional
%                  multiwavelet transform
%   thresh_decor2D scalar thresholding with decorrelation of 2-dimensional
%                  multiwavelet transform
%   thresh_vec2D   vector thresholding with decorrelation of 2-dimensional
%                  multiwavelet transform
%
%   covm1D_rrpe    computes covariance matrix of given 1-dimensional 
%                  multiwavelet transform with oversampled preprocessing
%   covm1D_appe    computes covariance matrix of given 1-dimensional 
%                  multiwavelet transform with critically sampled preprocessing
%   transop        computes transition operator of given multiwavelet transform(计算多重小波函数的转移矩阵)
%   multiplot      plots scaling functions and wavelets corresponding to the
%                  given multifilter bank(多重尺度函数和多重小波函数画图程序)
%
%   example1D      example how to use the package
%   simpcomp2D     simple compression of a 2-dimensional signal by retaining
%                  given number of largest coefficients
%   th2D           denoising of a 2-dimensional signal by various types of
%                  multiwavelet thresholding 
%
% To get more information about a _function_ type
%
%    help _function_
%
% List of related literature is in MWMP.README

