function [Pj,mu,Sigma,nu,Z,U,M]=tEMCore(Pj,mu,Sigma,nu,N,max_iteration,method, Input)
% EM steps for algorithms described in SS paper
% Sigma is cell array for unconstrained case, matrix for constrained case
% vector for diagonal case, cell array of vectors for idependent diagonal 
% covariances and empty for identity covariance.

%Initialization
format compact;
[NumRow NumCol] = size(Input);
%plot
% cla;
Loading = eye(NumCol);
% plot(Input(:,1), Input(:,2), 'k.', 'MarkerSize', 4);
% ph1 = plot(0,0);
% colors = repmat('b',20,1);
% hold on


g=length(Pj);
if g>1
   Pj=sparse(Pj);
end
rep=reshape(repmat(1:g,NumRow,1),g*NumRow,1);
rep_data=repmat(Input,g,1);
diffs=rep_data-mu(rep,:);%Distances to cluster center
delta=1;
Lmax=-inf;
L=Lmax;
nu_old=nu;
Lhist=[];
ghist=[];
SigmaNo=length(Sigma);
detSigma=sparse(zeros(1,SigmaNo));
for i=1:SigmaNo
   detSigma(i)=1/sqrt(det(Sigma{i}));      
end
M=zeros(NumRow,g);
for i=1:SigmaNo
   M(:,i)=sum(((sqrtm(pinv(Sigma{i}))*(diffs(find(rep==i),:)')).^2))'; %Mahalanobis distances         
end
c=gamma((nu+NumCol)/2)/(gamma(nu/2)*(pi*nu)^(NumCol/2));  
Prob=[c*exp(-(nu+NumCol)*log(1+M/nu)/2)*diag(detSigma)];%probabilities
k=0;%iteration counter
% Start Algorithm 
while (g>=1)&(k<max_iteration)
   while ((delta>0.1)|abs(nu-nu_old)>1e-1)&(k<max_iteration)  
      k=k+1;      
      %plot ellipses
%       delete(ph1)
%       ph1=[];
      
      for i=1:g
         sig1=Loading'*(Sigma{i}*Loading);%rotate covariance to Loading dimensions
         [V,D]=eig(sig1);
         center=(mu(i,:)-mean(Input))*Loading+mean(Input);
         if Pj(i)>0
%             ph1(i)=Ellipse(2*sqrt(D(1,1)),2*sqrt(D(2,2)),angle(V(1,1)+1i*V(2,1)),center(1),center(2),colors(i+2));
%             set(ph1(i),'LineWidth',2)%Pj(i)*20);
         else
%             ph1(i)=plot(center(1),center(2),'EraseMode','xor');
         end
      end
      drawnow   
      
      % E step      
      U=(nu+NumCol)./(nu+M);
      %temp=Prob*diag(Pj);
      %Z=diag(sparse(1./sum(temp,2)))*temp;
      Z=Prob*diag(Pj)./(sum(Prob*diag(Pj),2)*ones(1,g));
      ZU=Z.*U; 
      
      % M step
      switch method
      case 'agglomerate'
         deltaP=1;
         gtemp=g;
         while deltaP>10^-4
            Pjold=Pj;
            temp1=Prob*diag(Pj);
            temp2=sum(temp1,2);
            for j=1:g
               if Pj(j)>0
                  Pj(j)=max((sum(temp1(:,j)./temp2)-N/2)/(NumRow-N/2*gtemp),0);
                  if Pj(j)==0
                     temp2=sum(Prob*diag(Pj),2);
                     gtemp=gtemp-1;%number of components
                  end
               end
            end
            deltaP=norm(Pj-Pjold,1);
         end
      case 'regular'
         Pj=sparse(sum(Prob*diag(Pj)./(sum(Prob*diag(Pj),2)*ones(1,g)))/NumRow);
      end
      mu=(ZU'*Input)./(sum(ZU)'*ones(1,NumCol));
      %Update DOF parameter
      y=-sum(sum(Z.*(digamma((nu+NumCol)/2)+log(2./(nu+M))-U)))/sum(sum(Z));
      temp=1/(y+log(y)-1);
      nu_old=nu;
      erfInput = 0.6594*log(2.1971*temp);
      erfInput = real(erfInput); %MGChen
      nu=min(2*temp+0.0416*(1+erf(erfInput)),100);
      nu = real(nu); %MGChen
      %Calculate Covariance
      detSigma=[];
      for i=1:g
         diffs=Input-ones(NumRow,1)*mu(i,:);
         Sigma{i}=(((ZU(:,i)*ones(1,NumCol)).*diffs)'*diffs)/sum(ZU(:,i));
         detSigma(i)=sparse(1/sqrt(det(Sigma{i})));
         Sigma{i}(isnan(Sigma{i})) = 0;  %MGChen
         Sigma{i}(isinf(Sigma{i})) = 0;  %MGChen
         M(:,i)=sum(((sqrtm(pinv(Sigma{i}))*diffs').^2))'; %Mahalanobis distances 
      end
      % update Probabilities
      c=gamma((nu+NumCol)/2)/(gamma(nu/2)*(pi*nu)^(NumCol/2));  
      Prob=[c*exp(-(nu+NumCol)*log(1+M/nu)/2)*diag(detSigma)];
      switch method
      case 'agglomerate'
         % cycle is done, purge empty components
         tokeep=find(Pj);
         if length(tokeep)<g
            g=length(tokeep);
            rep=reshape(repmat(1:g,NumRow,1),g*NumRow,1);
            rep_data=repmat(Input,g,1);
            Pj=Pj(tokeep);
            mu=mu(tokeep,:);
            Prob=Prob(:,tokeep);
            Z=Z(:,tokeep);
            U=U(:,tokeep);
            ZU=ZU(:,tokeep);
            M=M(:,tokeep);
            Sigma=Sigma([tokeep]);
         end
         % Update and compare likelihood
         oldL=L;
         L=sum(log(sum(Prob*diag(sparse(Pj)),2)))-[N/2*sum(log(NumRow*Pj/12))+g/2*log(NumRow/12)+g*(N+1)/2]; %Mario expression 
         Lhist=[Lhist L];
         ghist=[ghist g];
         delta=abs(L-oldL); 
      end
   end
   switch method
   case 'agglomerate'
      if L>Lmax      
         Lmax=L;%current likelihood is optimal
         Pjopt=Pj;%store optimal parameters
         muoption=mu;
         Sigma_optimal=Sigma;
         nuopt=nu;
         Zopt=Z;
         Uopt=U;
         Mopt=M;
      else
         break
      end
      
      %purge smallest component
      [Y,ind]=min(Pj); %find minimal mixing proportion
      temp_I=1:g;
      I=setdiff(temp_I,ind);
      Pj=Pj(I);
      mu=mu(I,:);
      Prob=Prob(:,I);
      Z=Z(:,I);
      U=U(:,I);
      ZU=ZU(:,I);
      M=M(:,I);
      Sigma=Sigma([I]);
      g=g-1;   
      rep=reshape(repmat(1:g,NumRow,1),g*NumRow,1);
      rep_data=repmat(Input,g,1);
      delta=1;
   end
end

switch method
case 'agglomerate'
    if exist('muoption', 'var')==1
       mu=muoption;
       Sigma=Sigma_optimal;
       Pj=Pjopt;
       nu=nuopt;
       Z=Zopt;
       U=Uopt;
       M=Mopt;        
    end
end
mu=real(mu);
% Sigma=real(Sigma);
Pj=real(Pj);
nu=real(nu);
Z=real(Z);
U=real(U);
M=real(M); 