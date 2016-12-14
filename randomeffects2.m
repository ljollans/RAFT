function [ DATA Y braincorr outcomecorr]=randomeffects2(F, N, nClusters, MULTI1, MULTI2, MULTI3, effectsizes, frequencies, corrfreqs2use_clusters, corrs2use_clusters, corrfreqs2use_wholebrain, corrs2use_wholebrain);
% creating simulated brain data - Jollans et al., 2016
% %
% for 278 simulated data use this:
% [DATA Y braincorr outcomecorr]=randomeffects2(1000, 2000, 30, 3, .7, 1.5, [1], [100], [70 30], [0.9, 0.4], [100], [0.25]);
% 
% changing the hardcoded figures in this script will result in small
% changes to the correlation strengths between features and between
% features and the outcome variable. The numbers used here were arrived at
% by attempting to match the strentgh of feature-outcome correlations and
% feature-feature correlations as closely as possible to real neuroimaging
% data.
%
% this script uses the matlab built-in function mvnrnd as well as a number
% of other matlab functions
%
% for comments and questions please contact lee.jollans@gmail.com
%
% example usage: [X Y braincorr outcomecorr]=randomeffects(1000, 2000, 25);

basedata=rand(N, F);


b=shuffle(reversehist((frequencies*F)/100, effectsizes));
Y=basedata*b'; %we've now created a variable Y for which basedata explains 100% of the variance

%here we are creating an additional 'layer' with 'clusters' that have
%strong correlations within each cluster
clustersizes=sum_round(ones(nClusters,1)*F/nClusters);
ClusterM=[];
for cluster=1:nClusters
    %for each cluster we create a correlation matrix
    corr_frequencies=corrfreqs2use_clusters;
    corr_frequencies=(corr_frequencies*[(clustersizes(cluster)*(clustersizes(cluster)-1))/2])/100;
    corr_values=corrs2use_clusters;
    all_bb_corr=shuffle(reversehist(corr_frequencies, corr_values));
    Corr=ones(clustersizes(cluster), clustersizes(cluster));
    %fill in the correlation matrix for each variable in the cluster
    for n=1:clustersizes(cluster)
        Corr(n+1:end,n)=all_bb_corr(1:clustersizes(cluster)-n);
        Corr(n,n+1:end)=all_bb_corr(1:clustersizes(cluster)-n);
        all_bb_corr=all_bb_corr(clustersizes(cluster)-n+1:end);
    end
    % get the nearest Symmetric Positive Definite matrix to the correlation
    % matrix Corr to use as the covariance matrix for creating the
    % datapoints
    D=diag([ones(1,clustersizes(cluster))]);
    S=D*Corr*D;
    M=nearestSPD(S);
    data=mvnrnd(ones(1,clustersizes(cluster)), M, N);
    %Fill in the cluster layer
    ClusterM=[ClusterM, data(:,1:clustersizes(cluster))];
end

% next we make a layer of 'whole-brain' weakly correlated noise
corr_frequencies=corrfreqs2use_wholebrain;
corr_frequencies=(corr_frequencies*((F*(F-1))/2))/100;
corr_values=corrs2use_wholebrain;
all_bb_corr=shuffle(reversehist(corr_frequencies, corr_values));
Corr=ones(F, F);
%fill in the correlation matrix
for n=1:F
    Corr(n+1:end,n)=all_bb_corr(1:F-n);
    Corr(n,n+1:end)=all_bb_corr(1:F-n);
    all_bb_corr=all_bb_corr(F-n+1:end);
end
D=diag([ones(1,F)]);
S=D*Corr*D;
M=nearestSPD(S);
wholebrain_data=mvnrnd(ones(1,F), M, N);

% assemble everything
DATA=MULTI1*basedata+MULTI2*wholebrain_data+MULTI3*ClusterM;

% get the correlation strength
x=corr([DATA, Y]);
outcomecorr=abs(x(1:F,F+1)); %correlations between features and outcome
braincorr=x(1:F, 1:F); %correlations between features
braincorr(braincorr==1)=NaN;
braincorr=abs(unique(braincorr));
end

function A=reversehist(nelements, values)
total_elements=sum_round(nelements);
A=[];
for n=1:length(total_elements)
    A=[A,ones(1,total_elements(n))*values(n)];
end
end

function A=sum_round(B)
A=floor(B);
[val idx]=sort(B-A); val=fliplr(val); idx=fliplr(idx);
n=0;
while sum(A)<round(sum(B))
    n=n+1;
    A(idx(n))=A(idx(n))+1;
end
end

function X = shuffle(A)  
tmp=randperm(length(A));
X=A(tmp);
end

function Ahat = nearestSPD(A)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.

if nargin ~= 1
  error('Exactly one argument must be provided.')
end

% test for a square matrix A
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  % A was scalar and non-positive, so just return eps
  Ahat = eps;
  return
end

% symmetrize A into B
B = (A + A')/2;

% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[U,Sigma,V] = svd(B);
H = V*Sigma*V';

% get Ahat in the above formula
Ahat = (B+H)/2;

% ensure symmetry
Ahat = (Ahat + Ahat')/2;

% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
p = 1;
k = 0;
while p ~= 0
  [R,p] = chol(Ahat);
  k = k + 1;
  if p ~= 0
    % Ahat failed the chol test. It must have been just a hair off,
    % due to floating point trash, so it is simplest now just to
    % tweak by adding a tiny multiple of an identity matrix.
    mineig = min(eig(Ahat));
    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
end
end