function [L,H]=coef(flt)
%[L,H]=coef(flt)
%
%  This function returns coefficients of given scalar or matrix filter bank.
%  Notice that sometimes zero coefficients are added in the beginning of 
%  a filter. This ensures correct shift of the synthesis filter against 
%  analysis one. For details see [SW].
%                                                               
%  Input:                                                        
%    flt        string of characters, name of the filter; for admissible 
%               names and brief descriptions of filters see below 
%
%  Output:                                                       
%    L          这是低通滤波器的系数:
%               r by r*l real array, low-pass (scaling) filter;
%               r is the number of scaling functions,           
%               l is the number of terms in the dilation equation;
%               coefficients in L are organized as follows: L=[C1 C2 ... Cl]//C1...Cl是r*r的矩阵
%
%    H         这是高通滤波器的系数:
%               r by r*m real array, high-pass (wavelet) filter;
%               m is the number of terms in the wavelet equation;
%               coefficients in H are organized as follows: H=[D1 D2 ... Dm]
%
%  Admissible Names of the Filters and Brief Description of Correspoding Basis:
%   r is the number of scaling an wavelet functions,//r是维数
%   l is the number of scaling coefficients,        //l是低通滤波器系数的个数                  
%   m is the number of wavelet coefficients,        //m是高通滤波器系数的个数                  
%   A is the approximation order, see [HSS],        //A是逼近阶                  
%   S is the lower bound of Sobolev smoothness, see [D], [J])        
%
%   'haar'   Haar 2 coefficient orthogonal symmetric scalar filter bank 
%            (see [D]); r=1, l=2, m=2, A=1, S=0.4999
%   'd4'     Daubechies 4 coefficient orthogonal scalar filter bank (see [D]);
%            r=1, l=4, m=4, A=2, S=0.9999       
%   'la8'    Daubechies 8 coefficient least asymmetric orthogonal scalar 
%            filter bank (see [D]); r=1, l=8, m=8, A=4, S=1.7757
%
%%%%%%%%%以下是双正交小波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   'bi9'    9/7 coefficient symmetric biorthogonal  scalar filter bank,//这个小波常用来作为编码的比较标准 
%            dual to 'bi7' (see [D]); r=1, l=9, m=7, A=4, S=1.4101      
%   'bi7'    7/9 coefficient symmetric biorthogonal scalar filter bank, 
%            dual to 'bi9' (see [D]); r=1, l=7, m=9, A=4, S=2.1226
%
%
%   'bi5'    5/3 coefficient symmetric biorthogonal scalar filter bank, 
%            dual to 'bi3' (see [D]); r=1, l=5, m=3, A=2, S=0.4408
%   'bi3'    3/5 coefficient symmetric biorthogonal scalar filter bank, 
%            dual to 'bi5' (see [D]); r=1, l=3, m=5, A=2, S=1.4999
%
%%%%%%%%%以下是多重小波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   'ghm'    Geronimo-Hardin-Massopust orthogonal symmetric multi-filter bank 
%            (see [GHM], [SS]); r=2, l=4, m=4, A=2, S=1.4999
%   'cl'     Chui-Lian orthogonal symmetric  multi-filter bank (see [CL]);//这个多重小波长用来做编码 
%            r=2, l=3, m=3, A=2, S=1.0545 
%   'sa4'    orthogonal symmetric  multi-filter bank constructed by 
%            Shen, Tan, and Tham (see [STT]); r=2, l=4, m=4, A=1, S=0.9920
%
%%%%%%%以下是双正交多重小波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   'bih52s' most smooth biorthogonal symmetric multi-filter bank, 
%            dual to Hermite cubic filter bank 'bih32s' (see [TS]); 
%            r=2, l=5, m=3, A=2, S=0.8280
%   'bih32s' Hermite cubic multi-filter bank biorthogonal to 'bih52s' 
%            (see [TS]); r=2, l=3, m=5, A=4, S=2.4999 
%
%
%   'bih54n' non L2 biorthogonal symmetric multi-filter bank,
%            dual to Hermite cubic filter bank 'bih34n' (see [S]);
%            r=2, l=5, m=3, A=4, S=-0.6050
%   'bih34n' Hermite cubic multi-filter bank biorthogonal to  'bih54n' 
%            (see [S]); r=2, l=3, m=5, A=4, S=2.4999     
%
%
%   'bighm2' biorthogonal multi-filter bank obtained from 'ghm' filter bank 
%            by factoring out one approximation order, dual to 'bighm6' 
%            (see [S]); r=2, l=2, m=6, A=1, S=0.4999  %% This is for reconstruction 
%   'bighm6' biorthogonal multi-filter bank obtained from 'ghm' filter bank 
%            by adding one approximation order, dual to 'bighm2' (see [S]); 
%            r=2, l=6, m=2, A=3, S=2.4999  This is for decomosition 
%
%%%%%%以下是几个平衡多小波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   'cardbal2' orthogonal cardinal 2-balanced multi-filter bank constructed
%            by I. Selesnick (see [Se]); r=2, l=6, m=6, A=2, S=1.5261 
%   'cardbal3' orthogonal cardinal 3-balanced multi-filter bank constructed
%            by I. Selesnick (see [Se]); r=2, l=8, m=8, A=3, S=1.3345 
%   'cardbal4' orthogonal cardinal 4-balanced multi-filter bank constructed
%            by I. Selesnick (see [Se]); r=2, l=12, m=12, A=4, S=1.7979 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            以下是几个平衡多重小波：                        %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   'clbal'   this is balnance of Chui-Lian orthogonal symmetric  multi-filter bank;
%
%   'ghmbal'  this is balance of GHM multiwavelet
%
%
%
%  Example of Usage:
%   [L,H]=coef('bih52s')
 

% Author: Vasily Strela
% COPYRIGHT 1997,98 by Vasily Strela

if strcmp(flt,'haar')
  L=[1 1]/sqrt(2);
  H=[1 -1]/sqrt(2);

elseif strcmp(flt,'d4')
  L=[(1+sqrt(3))/(4*sqrt(2)) (3+sqrt(3))/(4*sqrt(2)) (3-sqrt(3))/(4*sqrt(2)) (1-sqrt(3))/(4*sqrt(2))];
  H=[(1-sqrt(3))/(4*sqrt(2)) (-3+sqrt(3))/(4*sqrt(2)) (3+sqrt(3))/(4*sqrt(2)) (-1-sqrt(3))/(4*sqrt(2))];

elseif strcmp(flt,'la8')
  L=[-0.07576571478950 -0.02963552764600   0.49761866763278   0.80373875180513 0.29785779560530  -0.09921954357663  -0.01260396726203   0.03222310060405];
  H=[-0.03222310060405   -0.01260396726203  0.09921954357663  0.29785779560530 -0.80373875180513  0.49761866763278  0.02963552764600   -0.07576571478950];

elseif strcmp(flt,'bi9')
  L=[0.03782845550700 -0.02384946501938  -0.11062440441842  0.37740285561265 0.85269867900940  0.37740285561265  -0.11062440441842  -0.02384946501938 0.03782845550700];
  H=[0 0 0.06453888262894  -0.04068941760956 -0.41809227322221 0.78848561640566 -0.41809227322221  -0.04068941760956  0.06453888262894];

elseif strcmp(flt,'bi7')
  L=[0 -0.06453888262894  -0.04068941760956   0.41809227322221   0.78848561640566 0.41809227322221  -0.04068941760956  -0.06453888262894];
  H=[0  0.03782845550700   0.02384946501938  -0.11062440441842  -0.37740285561265 0.85269867900940  -0.37740285561265  -0.11062440441842  0.02384946501938 0.03782845550700]; 

elseif strcmp(flt,'bi5')
  L=[-1/4 2/4 6/4 2/4 -1/4]/sqrt(2);
  H=[0 0 -1/4 2/4 -1/4]/sqrt(2);

elseif strcmp(flt,'bi3')
  L=[0 1/4 2/4 1/4]*sqrt(2);
  H=[0 -1/4 -2/4 6/4 -2/4 -1/4]*sqrt(2);

elseif strcmp(flt,'ghm')
  L=[ 3/(5*sqrt(2))  4/5            3/(5*sqrt(2)) 0          0     0               0    0
     -1/20          -3/(10*sqrt(2)) 9/20          1/sqrt(2)  9/20 -3/(10*sqrt(2)) -1/20 0];

  H=[-1/20           -3/(10*sqrt(2))  9/20           -1/sqrt(2) 9/20             -3/(10*sqrt(2)) -1/20           0
      1/(10*sqrt(2))  3/10           -9/(10*sqrt(2))  0         9/(10*sqrt(2))   -3/10           -1/(10*sqrt(2)) 0];

elseif strcmp(flt,'cl')
  L=[1/(2*sqrt(2))       -1/(2*sqrt(2))       1/sqrt(2) 0              1/(2*sqrt(2))        1/(2*sqrt(2))
     sqrt(7)/(4*sqrt(2)) -sqrt(7)/(4*sqrt(2)) 0         1/(2*sqrt(2)) -sqrt(7)/(4*sqrt(2)) -sqrt(7)/(4*sqrt(2))];

  H=[1/(2*sqrt(2)) -1/(2*sqrt(2)) -1/sqrt(2) 0                   1/(2*sqrt(2)) 1/(2*sqrt(2))
   -1/(4*sqrt(2))  1/(4*sqrt(2))  0         sqrt(7)/(2*sqrt(2)) 1/(4*sqrt(2)) 1/(4*sqrt(2))];

elseif strcmp(flt,'bih52s')
  L=[-73/648  -77/972     1/2     89/486 397/324 0         1/2    -89/486 -73/648 77/972
      773/1080 3229/6480 -187/60 -91/81   0      6091/3240 187/60 -91/81  -773/1080 3229/6480]/sqrt(2);
  H=[0 0 0 0 -1/4 -1/8 1/2 0   -1/4 1/8
     0 0 0 0  3/8  1/8 0   1/2 -3/8 1/8]/sqrt(2);

elseif strcmp(flt,'bih32s')
  L=[0 0 1/2  3/4 1  0   1/2 -3/4
     0 0 -1/8 -1/8 0  1/2 1/8 -1/8]/sqrt(2);
  H=[0 0 67/240  7/240 -1     -187/60 173/120 0    -1      187/60 67/240 -7/240
     0 0 -95/972 -1/162  89/243 91/81  0       26/9 -89/243 91/81  95/972 -1/162]/sqrt(2);

elseif strcmp(flt,'bih54n')
  L=[-3/16 -5/48  1/2   1/12 11/8  0   1/2  -1/12 -3/16  5/48
     33/32  9/16 -21/8 -3/8   0    5/2 21/8 -3/8  -33/32 9/16]/sqrt(2);
  H=[0 0 0 0  1/2  1/4 -1  0   1/2 -1/4
     0 0 0 0 -3/8 -1/8  0 -1/2 3/8 -1/8]/sqrt(2);

elseif strcmp(flt,'bih34n')
  L=[0 0 1/2  3/4 1  0   1/2 -3/4
     0 0 -1/8 -1/8 0  1/2 1/8 -1/8]/sqrt(2);
  H=[0 0 -5/64  3/64  1/2  21/16 -27/32 0     1/2 -21/16 -5/64 -3/64
     0 0  1/96 -1/32 -1/6 -3/8    0    -57/16 1/6 -3/8   -1/96 -1/32]/sqrt(2);

elseif strcmp(flt,'sa4')
  L=[(32 + 8*sqrt(15))^(-1),  1/8,  (31 + 8*sqrt(15))/(8*(4 + sqrt(15))), 1/8, (31 + 8*sqrt(15))/(8*(4 + sqrt(15))), -1/8, (32 + 8*sqrt(15))^(-1), -1/8;
     (32 + 8*sqrt(15))^(-1), -1/8, -(31 + 8*sqrt(15))/(8*(4 + sqrt(15))), 1/8, (31 + 8*sqrt(15))/(8*(4 + sqrt(15))),  1/8, -1/(8*(4 + sqrt(15))),  -1/8]/sqrt(2);

  H=[-1/8, (32 + 8*sqrt(15))^(-1), 1/8, -(31 + 8*sqrt(15))/(8*(4 + sqrt(15))), 1/8, (31 + 8*sqrt(15))/(8*(4 + sqrt(15))),  -1/8, -1/(8*(4 + sqrt(15)));
     -1/8, -1/(8*(4 + sqrt(15))), -1/8, -(31 + 8*sqrt(15))/(8*(4 + sqrt(15))), 1/8, -(31 + 8*sqrt(15))/(8*(4 + sqrt(15))), 1/8, -1/(8*(4 + sqrt(15)))]/sqrt(2);

elseif strcmp(flt,'bighm2')
  L=[0,0,0,0,1,   0,   1,   0;
     0,0,0,0,8/5,-2/5,-8/5,-2/5]/sqrt(2);

  H=[3/13,-2/13,-1,    -2/13,10/13,40/13,10/13,-40/13,-1,    2/13,3/13,2/13;   
     3/2, -1,   -13/2, -1,  -4,    26,   4,     26,    13/2,-1,  -3/2,-1]/sqrt(2);

elseif strcmp(flt,'bighm6')
  L=[-1/40,1/40,1/40,-9/40,1,   -1/4,  1,    1/4,  1/40, 9/40,-1/40,-1/40;
     -1/40,1/40,1/40,-9/40,13/20,1/10,-13/20,1/10,-1/40,-9/40, 1/40, 1/40]/sqrt(2);

  H=[0,0,0,0,0,    13/40,0,   -13/40;
     0,0,0,0,1/100,1/25,-1/100,1/25]/sqrt(2);

elseif strcmp(flt,'cardbal2')
L=[0.02209708691208                 0   0.17396999725850   0.70710678118655 ...
0.66291260736239                  0  -0.17116329922036                  0 ...
0.02209708691208                  0  -0.00280669803814                  0;
0.00280669803814                  0   0.02209708691208                  0 ...
0.17116329922036   0.70710678118655   0.66291260736239                  0 ...
-0.17396999725850                  0   0.02209708691208                  0];

H=[-0.02209708691208                0  -0.17396999725850   0.70710678118655 ...
-0.66291260736239                  0   0.17116329922036                  0 ...
-0.02209708691208                  0   0.00280669803814                  0;
-0.00280669803814                  0  -0.02209708691208                  0 ...
-0.17116329922036   0.70710678118655  -0.66291260736239                  0 ...
0.17396999725850                  0  -0.02209708691208                  0];

elseif strcmp(flt,'cardbal3')
L=[0.00777141529460                  0   0.01890191323078                 0 ...
0.15346244941283   0.70710678118655   0.67249812840628                  0 ...
-0.15346244941283                  0   0.01251156586818                  0 ...
-0.00777141529460                  0   0.00319517368130                  0;
0.00319517368130                  0   0.00777141529460                  0 ...
0.01251156586818                  0   0.15346244941283   0.70710678118655 ...
0.67249812840628                  0  -0.15346244941283                  0 ...
0.01890191323078                  0  -0.00777141529460                  0];

H=[-0.00777141529460                  0  -0.01890191323078                0 ...
-0.15346244941283   0.70710678118655  -0.67249812840628                  0 ...
 0.15346244941283                  0  -0.01251156586818                  0 ...
0.00777141529460                  0  -0.00319517368130                  0;
-0.00319517368130                  0  -0.00777141529460                  0 ...
-0.01251156586818                  0  -0.15346244941283   0.70710678118655 ...
-0.67249812840628                  0   0.15346244941283                  0 ...
-0.01890191323078                  0   0.00777141529460                  0];

elseif strcmp(flt,'cardbal4')
L=[-0.00000045995944                  0  -0.00172617386509               0 ...
-0.00491803042978                  0   0.02900162132252                  0 ...
0.19092914488214   0.70710678118655   0.65255620837152                  0 ...
-0.18972126188006                  0   0.02900081607293                  0 ...
0.00310620592666                  0  -0.00172552966542                  0 ...
0.00060440146048                  0  -0.00000016104992                  0;
-0.00000016104992                  0  -0.00060440146048                  0 ...
-0.00172552966542                  0  -0.00310620592666                  0 ...
0.02900081607293                  0   0.18972126188006   0.70710678118655 ...
0.65255620837152                  0  -0.19092914488214                  0 ...
0.02900162132252                  0   0.00491803042978                  0 ...
-0.00172617386509                  0   0.00000045995944                  0];

H=[0.00000045995944                  0   0.00172617386509                0 ...
0.00491803042978                  0  -0.02900162132252                  0 ...
-0.19092914488214   0.70710678118655  -0.65255620837152                  0 ...
0.18972126188006                  0  -0.02900081607293                  0 ...
-0.00310620592666                  0   0.00172552966542                  0 ...
-0.00060440146048                  0   0.00000016104992                  0;
0.00000016104992                  0   0.00060440146048                  0 ...
0.00172552966542                  0   0.00310620592666                  0 ...
-0.02900081607293                  0  -0.18972126188006   0.70710678118655 ...
-0.65255620837152                  0   0.19092914488214                  0 ...
-0.02900162132252                  0  -0.00491803042978                  0 ...
0.00172617386509                  0  -0.00000045995944                  0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   以下是平衡多重小波滤波器         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(flt,'clbal')   %%%%     这是cl平衡多重小波，这常用来编码    %%%%%
R=[1 -1;1 1];

C11=[1/(2*sqrt(2))  -1/(2*sqrt(2)) ; sqrt(7)/(4*sqrt(2))   -sqrt(7)/(4*sqrt(2))];
C21=[1/sqrt(2)      0              ; 0                      1/(2*sqrt(2))      ];
C31=[1/(2*sqrt(2))  1/(2*sqrt(2))  ; -sqrt(7)/(4*sqrt(2))  -sqrt(7)/(4*sqrt(2))];

D11=[1/(2*sqrt(2)) -1/(2*sqrt(2))  ; -1/(4*sqrt(2))        1/(4*sqrt(2))       ];
D21=[-1/sqrt(2)    0               ; 0                     sqrt(7)/(2*sqrt(2)) ];
D31=[1/(2*sqrt(2)) 1/(2*sqrt(2))   ; 1/(4*sqrt(2))         1/(4*sqrt(2))       ];

C1=R*C11*R';C2=R*C21*R';C3=R*C31*R';
D1=R*D11*R';D2=R*D21*R';D3=R*D31*R';

L=[C1,C2,C3];
H=[D1,D2,D3];

elseif strcmp(flt,'ghmbal')  %%%%%%%%%%   这是ghm平衡多重小波   %%%%%%%%%%%%%%
R=[1 -1;1 1];

C11=[3/(5*sqrt(2))  4/5;-1/20   -3/(10*sqrt(2))];
C21=[3/(5*sqrt(2))  0  ;9/20          1/sqrt(2)];
C31=[0              0  ;9/20    -3/(10*sqrt(2))];
C41=[0              0  ; -1/20  0              ];

C1=R*C11*R';C2=R*C21*R';C3=R*C31*R';C4=R*C41*R';

D11=[-1/20           -3/(10*sqrt(2));1/(10*sqrt(2))  3/10 ];
D21=[9/20           -1/sqrt(2)      ;-9/(10*sqrt(2))  0   ];
D31=[9/20            -3/(10*sqrt(2));9/(10*sqrt(2))   -3/10];
D41=[-1/20           0              ;-1/(10*sqrt(2)) 0];

D1=R*D11*R';D2=R*D21*R';D3=R*D31*R';D4=R*D41*R';

L=[C1,C2,C3,C4];
H=[D1,D2,D3,D4];

elseif strcmp(pflt,'opt')         %% The multiwavelets is optiminal_reconstruction 
    R=(1/sqrt(2))*[1 -1 ;1 1];    %% It seems some problem
    
    C11=0.5*[1 0 ;-sqrt(3)/2 1/2];
    C21=0.5*[1 0 ; sqrt(3)/2 1/2];
    
    C1=R*C11*R';C2=R*C21*R';
    
    D11=0.5*[-1/2 -sqrt(3)/2 ; 0  -1];
    D21=0.5*[ 1/2 -sqrt(3)/2 ; 0   1];
    
    D1=R*D11*R';D2=R*D21*R';
    
    L=[C1,C2];
    H=[D1,D2];
    

end 

