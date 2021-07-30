function [Class Model] = tdistEM(Input, PenaltyFactor)
% t-distribution EM proposed by Shy Shoham, J Neurosci Methods2003
% Adapted from SAC by Shy Shoham
% [~, Input] = princomp(Input);
% Input = Input(:,1:3);
% Initialization
Method = 't';
[NumRow NumCol] = size(Input);
Loading = eye(NumCol);
NumCluster = NumCol*15;
NumCluster=NumCluster/PenaltyFactor;
NumFuzzyCluster = 10;
while round(NumRow/NumCluster)<NumFuzzyCluster
    NumCluster = NumCluster/2;
    NumFuzzyCluster = NumFuzzyCluster-2;
end
[Center FuzzyMatrix] = fcm(Input, NumFuzzyCluster, [2,20,1,0]);
FuzzyMatrix=FuzzyMatrix';
% Estimate error covariance matrix
Rep = reshape(repmat(1:NumFuzzyCluster, NumRow, 1), NumFuzzyCluster*NumRow, 1);
Dist = repmat(Input, NumFuzzyCluster, 1)-Center(Rep,:);
for i = 1:NumFuzzyCluster
   SD{i} = (((FuzzyMatrix(:,i)*ones(1, NumCol)).*Dist(find(Rep==i),:))'*Dist(find(Rep==i),:))/sum(FuzzyMatrix(:,i));
end
Pj = [(sum(FuzzyMatrix))/sum(sum(FuzzyMatrix))];
% EM clustering
NU=20;  %initial nu value
[Pj,Center,SD,NU,Z,FuzzyMatrix,M] = ...
    tEMCore(Pj,Center,SD,NU,NumCluster,10,'regular',Input);  %10 iterations of regular EM
[Pj,Center,SD,NU,Z,FuzzyMatrix,M] = ...
    tEMCore(Pj,Center,SD,NU,NumCluster,500,'agglomerate',Input);  %modified Mario algorithm
[Pj,Center,SD,NU,Z,FuzzyMatrix,M] = ...
    tEMCore(Pj,Center,SD,NU,NumCluster,10,'regular',Input);  %10 iterations of regular EM
NumFuzzyCluster = length(Pj);
ZU = Z.*FuzzyMatrix;
% Max_M = InverseDistribution(NU, NumCol, 0.998);
x = [0.1:.2:NumCol*5];
y=cumsum(betad(1./(1+x/NU),NU/2+2,NumCol/2));
y=y/y(end);
Max_M = spline(y,x,0.998);

Outliers = find(min(M, 2)>Max_M);
Center = Center*(Loading(:,1:NumCol)')+ones(size(Center,1),1)*mean(Input);
for i = 1:NumFuzzyCluster
  Dist = Input-ones(NumRow,1)*Center(i,:);
  SD{i} = (((ZU(:,i)*ones(1,NumCol)).*Dist)'*Dist)/sum(ZU(:,i));
end
ClusterInfo.Center = Center;
ClusterInfo.Sigma = SD;
ClusterInfo.Proportion = Pj;
[Y Class] = max((Z)', [], 1);
ClusterInfo.Units=[];
if NumFuzzyCluster>1 %decide which clusters contain local field potentials or overlaps
   for i = 1:NumFuzzyCluster
      sig1 = Loading(:,1:2)'*(ClusterInfo.Sigma{i}*Loading(:,1:2));
      f(i) = prod(diag(sig1));
   end
   Garbage = find(f>10*median(f)); %large covariance indicates a garbage collector
   NonGarbage = setdiff(1:NumFuzzyCluster, Garbage);
   [~, NewOrder] = sort(sum(detrend(ClusterInfo.Center(NonGarbage,:)').^2));  %sort by ascending energy
   sig1 = Loading(:,1:2)'*(ClusterInfo.Sigma{NonGarbage(NewOrder(1))}*Loading(:,1:2));
   mu1 = ClusterInfo.Center(NonGarbage(NewOrder(1)),:)*Loading(:,1:2);
   if sum(mu1.^2)<4*sum(diag(sig1)) % LFP decision
      NumUnit = length(NewOrder)-1;
      ClusterInfo.Units(NonGarbage(NewOrder)) = 0:NumUnit;
   else
      NumUnit = length(NewOrder);
      ClusterInfo.Units(NonGarbage(NewOrder)) = 1:NumUnit;
   end
   ClusterInfo.Units(Garbage) = 255;
else
   NumUnit = 1;
   ClusterInfo.Units = 1;
end
Class = ClusterInfo.Units(Class)+1;
Class(Outliers) = 0;
if strcmpi(Method, 't')
   ClusterInfo.nu = NU;
else
   ClusterInfo.nu = inf;
end
[~, Index] = sort(ClusterInfo.Units);
Class = Class';
Model.Center = ClusterInfo.Center(Index,:);
Model.Sigma = ClusterInfo.Sigma(Index);