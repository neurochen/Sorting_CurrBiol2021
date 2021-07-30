% function dum=example1D(n0)
%dum=example1D(n0)
%
%  This function is an example of usage of all scalar and multiwavelet 
%  transforms built into MWMP. 
%  A random signal of given length is decomposed and then reconstructed.
%  Root mean square error between initial and reconstructed signals is computed.
%  For description of the names of the filters see coef.m and coef_prep.m 
%  
%  Input:                                                       
%    n0         integer, length of the signal;
%               n0 must be of the form constant*2^k, k > 3
%
%  rr 是重复行预滤波；
%  appe 是给定逼近阶预滤波；
%  id 这是用于平衡小波预滤波；
%
%      对于1D的信号的分解显示:先对输入信号进行预滤波,然后再进行多小波分解,对分解后的
%  频带信号可以进行后滤波,之后可以将分解后的频带显示．
%
%
%
%  Example of Usage:
%    example1D(384)

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela
n0 = 48;
maxlevel=round(log(n0)/log(2))-3;
f=randn(1,n0);
subplot(511); plot(f);

'Decomposition: ghm, Reconstruction: ghm, Preprocessing: rr'
fp=prep1D_rr(f,'ghm');
subplot(512); plot(fp');
transf=dec1D_pe(fp,'ghm',maxlevel);
subplot(513); plot(transf');
fhatp=rec1D_pe(transf,'ghm',maxlevel);
subplot(514); plot(fhatp');
fhat=postp1D_rr(fhatp,'ghm');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: ghm, Reconstruction: ghm, Preprocessing: ghmap'
fp=prep1D_appe(f,'ghmap');
transf=dec1D_pe(fp,'ghm',maxlevel-1);
fhatp=rec1D_pe(transf,'ghm',maxlevel-1);
fhat=postp1D_appe(fhatp,'ghmap');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: cl, Reconstruction: cl, Preprocessing: rr'
fp=prep1D_rr(f,'cl');
transf=dec1D_pe(fp,'cl',maxlevel);
fhatp=rec1D_pe(transf,'cl',maxlevel);
fhat=postp1D_rr(fhatp,'cl');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: cl, Reconstruction: cl, Preprocessing: clap'
fp=prep1D_appe(f,'clap');
transf=dec1D_pe(fp,'cl',maxlevel-1);
fhatp=rec1D_pe(transf,'cl',maxlevel-1);
fhat=postp1D_appe(fhatp,'clap');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: sa4, Reconstruction: sa4, Preprocessing: rr'
fp=prep1D_rr(f,'sa4');
transf=dec1D_pe(fp,'sa4',maxlevel-1);
fhatp=rec1D_pe(transf,'sa4',maxlevel-1);
fhat=postp1D_rr(fhatp,'sa4');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: sa4, Reconstruction: sa4, Preprocessing: sa4ap'
fp=prep1D_appe(f,'sa4ap');
transf=dec1D_pe(fp,'sa4',maxlevel-1);
fhatp=rec1D_pe(transf,'sa4',maxlevel-1);
fhat=postp1D_appe(fhatp,'sa4ap');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: cardbal2, Reconstruction: cardbal2, Preprocessing: id'
fp=prep1D_appe(f,'id');
transf=dec1D_pe(fp,'cardbal2',maxlevel-1);
fhatp=rec1D_pe(transf,'cardbal2',maxlevel-1);
fhat=postp1D_appe(fhatp,'id');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: cardbal4, Reconstruction: cardbal4, Preprocessing: id'

'Decomposition: cardbal3, Reconstruction: cardbal3, Preprocessing: id'
fp=prep1D_appe(f,'id');
transf=dec1D_pe(fp,'cardbal3',maxlevel-1);
fhatp=rec1D_pe(transf,'cardbal3',maxlevel-1);
fhat=postp1D_appe(fhatp,'id');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: cardbal4, Reconstruction: cardbal4, Preprocessing: id'
fp=prep1D_appe(f,'id');
transf=dec1D_pe(fp,'cardbal4',maxlevel-1);
fhatp=rec1D_pe(transf,'cardbal4',maxlevel-1);
fhat=postp1D_appe(fhatp,'id');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bih52s, Reconstruction: bih32s, Preprocessing: rr'
fp=prep1D_rr(f,'bih52s');
transf=dec1D_pe(fp,'bih52s',maxlevel);
fhatp=rec1D_pe(transf,'bih32s',maxlevel);
fhat=postp1D_rr(fhatp,'bih52s');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bih52s, Reconstruction: bih32s, Preprocessing: bih5ap'
fp=prep1D_appe(f,'bih5ap');
transf=dec1D_pe(fp,'bih52s',maxlevel-1);
fhatp=rec1D_pe(transf,'bih32s',maxlevel-1);
fhat=postp1D_appe(fhatp,'bih5ap');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bih54n, Reconstruction: bih34n, Preprocessing: rr'
fp=prep1D_rr(f,'bih54n');
transf=dec1D_pe(fp,'bih54n',maxlevel);
fhatp=rec1D_pe(transf,'bih34n',maxlevel);
fhat=postp1D_rr(fhatp,'bih54n');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bih54n, Reconstruction: bih34n, Preprocessing: bih5ap'
fp=prep1D_appe(f,'bih5ap');
transf=dec1D_pe(fp,'bih54n',maxlevel-1);
fhatp=rec1D_pe(transf,'bih34n',maxlevel-1);
fhat=postp1D_appe(fhatp,'bih5ap');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bighm2, Reconstruction: bighm6, Preprocessing: rr'
fp=prep1D_rr(f,'bighm2');
transf=dec1D_pe(fp,'bighm2',maxlevel);
fhatp=rec1D_pe(transf,'bighm6',maxlevel);
fhat=postp1D_rr(fhatp,'bighm2');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bighm2, Reconstruction: bighm6, Preprocessing: bighm2ap'
fp=prep1D_appe(f,'bighm2ap');
transf=dec1D_pe(fp,'bighm2',maxlevel-1);
fhatp=rec1D_pe(transf,'bighm6',maxlevel-1);
fhat=postp1D_appe(fhatp,'bighm2ap');
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: haar, Reconstruction: haar, Preprocessing: none'
transf=dec1D_pe(f,'haar',maxlevel);
fhat=rec1D_pe(transf,'haar',maxlevel);
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: d4, Reconstruction: d4, Preprocessing: none'
transf=dec1D_pe(f,'d4',maxlevel);
fhat=rec1D_pe(transf,'d4',maxlevel);
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: la8, Reconstruction: la8, Preprocessing: none'
transf=dec1D_pe(f,'la8',maxlevel);
fhat=rec1D_pe(transf,'la8',maxlevel);
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bi9, Reconstruction: bi7, Preprocessing: none'
transf=dec1D_pe(f,'bi9',maxlevel);
fhat=rec1D_pe(transf,'bi7',maxlevel);
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0

'Decomposition: bi5, Reconstruction: bi3, Preprocessing: none'
transf=dec1D_pe(f,'bi5',maxlevel);
fhat=rec1D_pe(transf,'bi3',maxlevel);
subplot(515); plot(fhat);
rmse_between_initial_and_reconstructed=norm(f-fhat,'fro')/n0