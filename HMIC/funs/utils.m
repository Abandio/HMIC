

function res = obj_value(X,A,Z,P,View,alpha,beta,delta)
res = 0;
sum_Z = 0;
    for i = 1:View
    for j = 1:View
        if i ~= j
            res = res - alpha * trace(Z{i}' * Z{i} * Z{j}' * Z{j});
        end
    end

    res = res + norm(X{i} - Z{i}*A{i},'fro')^2 ;
    sum_Z = sum_Z + delta(i) * Z{i};
    end
    res = res + beta * norm(sum_Z - P,'fro')^2;
end


%  min  x'*A*x - x'*b
%  s.t. x'1=1, x>=0
function [x, obj]=SimplexQP_acc(A, b, x0)


NIter = 500;
NStop = 20;

[n] = size(A,1);
if nargin < 3
    x = 1/n*ones(n,1);
else
    x = x0;
end

x1 = x; 
t = 1;
t1 = 0;
r = 0.5;  % r=1/mu;
%obj = zeros(NIter,1);
for iter = 1:NIter
    p = (t1-1)/t;
    s = x + p*(x-x1); 
    x1 = x;
    g = 2*A*s - b;
    ob1 = x'*A*x - x'*b;
    for it = 1:NStop
        z = s - r*g;
        z = EProjSimplex_new(z,1);  % z'1=1;z>=0; z=alpha;
        ob = z'*A*z - z'*b;
        if ob1 < ob
            r = 0.5*r; % rho=2;
        else
            break;
        end
    end
    if it == NStop
        obj(iter) = ob;
        %disp('not');
        break;
    end
    x = z;
    t1 = t;
    t = (1+sqrt(1+4*t^2))/2;
    
    
    obj(iter) = ob;
end
end



function Z = my_KNN(data,marks,opts)

if (~exist('opts','var'))
   opts = [];
end

r = 5;
if isfield(opts,'r')
    r = opts.r;
end


p = size(marks,1);
if isfield(opts,'p')
    p = opts.p;
end


nSmp=size(data,1);

% Z construction
D = EuDist2(data,marks,0);

if isfield(opts,'sigma')
    sigma = opts.sigma;
else
    sigma = mean(mean(D));
end

dump = zeros(nSmp,r);
idx = dump;
for i = 1:r
    [dump(:,i),idx(:,i)] = min(D,[],2);
    temp = (idx(:,i)-1)*nSmp+[1:nSmp]';
    D(temp) = 1e100;
end

dump = exp(-dump/(2*sigma^2));
%dump = dump./repmat(sum(dump,2),1,size(dump,2));

sumD = sum(dump,2);
Gsdx = bsxfun(@rdivide,dump,sumD);
Gidx = repmat([1:nSmp]',1,r);
Gjdx = idx;
Z=sparse(Gidx(:),Gjdx(:),Gsdx(:),nSmp,p);


end


function D = EuDist2(fea_a,fea_b,bSqrt)
%EUDIST2 Efficiently Compute the Euclidean Distance Matrix by Exploring the
%Matlab matrix operations.
%
%   D = EuDist(fea_a,fea_b)
%   fea_a:    nSample_a * nFeature
%   fea_b:    nSample_b * nFeature
%   D:      nSample_a * nSample_a
%       or  nSample_a * nSample_b
%
%    Examples:
%
%       a = rand(500,10);
%       b = rand(1000,10);
%
%       A = EuDist2(a); % A: 500*500
%       D = EuDist2(a,b); % D: 500*1000
%
%   version 2.1 --November/2011
%   version 2.0 --May/2009
%   version 1.0 --November/2005
%
%   Written by Deng Cai (dengcai AT gmail.com)


if ~exist('bSqrt','var')
    bSqrt = 1;
end

if (~exist('fea_b','var')) || isempty(fea_b)
    aa = sum(fea_a.*fea_a,2);
    ab = fea_a*fea_a';
    
    if issparse(aa)
        aa = full(aa);
    end
    
    D = bsxfun(@plus,aa,aa') - 2*ab;
    D(D<0) = 0;
    if bSqrt
        D = sqrt(D);
    end
    D = max(D,D');
else
    aa = sum(fea_a.*fea_a,2);
    bb = sum(fea_b.*fea_b,2);
    ab = fea_a*fea_b';

    if issparse(aa)
        aa = full(aa);
        bb = full(bb);
    end

    D = bsxfun(@plus,aa,bb') - 2*ab;
    D(D<0) = 0;
    if bSqrt
        D = sqrt(D);
    end
end


end


function [label, center, bCon, sumD, D] = litekmeans(X, k, varargin)
%LITEKMEANS K-means clustering, accelerated by matlab matrix operations.
%
%   label = LITEKMEANS(X, K) partitions the points in the N-by-P data matrix
%   X into K clusters.  This partition minimizes the sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  Rows of X correspond to points, columns correspond to
%   variables.  KMEANS returns an N-by-1 vector label containing the
%   cluster indices of each point.
%
%   [label, center] = LITEKMEANS(X, K) returns the K cluster centroid
%   locations in the K-by-P matrix center.
%
%   [label, center, bCon] = LITEKMEANS(X, K) returns the bool value bCon to
%   indicate whether the iteration is converged.  
%
%   [label, center, bCon, SUMD] = LITEKMEANS(X, K) returns the
%   within-cluster sums of point-to-centroid distances in the 1-by-K vector
%   sumD.    
%
%   [label, center, bCon, SUMD, D] = LITEKMEANS(X, K) returns
%   distances from each point to every centroid in the N-by-K matrix D. 
%
%   [ ... ] = LITEKMEANS(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by KMEANS.  Parameters are:
%
%   'Distance' - Distance measure, in P-dimensional space, that KMEANS
%      should minimize with respect to.  Choices are:
%            {'sqEuclidean'} - Squared Euclidean distance (the default)
%             'cosine'       - One minus the cosine of the included angle
%                              between points (treated as vectors). Each
%                              row of X SHOULD be normalized to unit. If
%                              the intial center matrix is provided, it
%                              SHOULD also be normalized.
%
%   'Start' - Method used to choose initial cluster centroid positions,
%      sometimes known as "seeds".  Choices are:
%         {'sample'}  - Select K observations from X at random (the default)
%          'cluster' - Perform preliminary clustering phase on random 10%
%                      subsample of X.  This preliminary phase is itself
%                      initialized using 'sample'. An additional parameter
%                      clusterMaxIter can be used to control the maximum
%                      number of iterations in each preliminary clustering
%                      problem.
%           matrix   - A K-by-P matrix of starting locations; or a K-by-1
%                      indicate vector indicating which K points in X
%                      should be used as the initial center.  In this case,
%                      you can pass in [] for K, and KMEANS infers K from
%                      the first dimension of the matrix.
%
%   'MaxIter'    - Maximum number of iterations allowed.  Default is 100.
%
%   'Replicates' - Number of times to repeat the clustering, each with a
%                  new set of initial centroids. Default is 1. If the
%                  initial centroids are provided, the replicate will be
%                  automatically set to be 1.
%
% 'clusterMaxIter' - Only useful when 'Start' is 'cluster'. Maximum number
%                    of iterations of the preliminary clustering phase.
%                    Default is 10.  
%
%
%    Examples:
%
%       fea = rand(500,10);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50);
%
%       fea = rand(500,10);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50, 'Replicates', 10);
%
%       fea = rand(500,10);
%       [label, center, bCon, sumD, D] = litekmeans(fea, 5, 'MaxIter', 50);
%       TSD = sum(sumD);
%
%       fea = rand(500,10);
%       initcenter = rand(5,10);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50, 'Start', initcenter);
%
%       fea = rand(500,10);
%       idx=randperm(500);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50, 'Start', idx(1:5));
%
%
%   See also KMEANS
%
%    [Cite] Deng Cai, "Litekmeans: the fastest matlab implementation of
%           kmeans," Available at:
%           http://www.zjucadcg.cn/dengcai/Data/Clustering.html, 2011. 
%
%   version 2.0 --December/2011
%   version 1.0 --November/2011
%
%   Written by Deng Cai (dengcai AT gmail.com)


if nargin < 2
    error('litekmeans:TooFewInputs','At least two input arguments required.');
end

[n, p] = size(X);


pnames = {   'distance' 'start'   'maxiter'  'replicates' 'onlinephase' 'clustermaxiter'};
dflts =  {'sqeuclidean' 'sample'       []        []        'off'              []        };
[eid,errmsg,distance,start,maxit,reps,online,clustermaxit] = getargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('litekmeans:%s',eid),errmsg);
end

if ischar(distance)
    distNames = {'sqeuclidean','cosine'};
    j = strcmpi(distance, distNames);
    j = find(j);
    if length(j) > 1
        error('litekmeans:AmbiguousDistance', ...
            'Ambiguous ''Distance'' parameter value:  %s.', distance);
    elseif isempty(j)
        error('litekmeans:UnknownDistance', ...
            'Unknown ''Distance'' parameter value:  %s.', distance);
    end
    distance = distNames{j};
else
    error('litekmeans:InvalidDistance', ...
        'The ''Distance'' parameter value must be a string.');
end


center = [];
if ischar(start)
    startNames = {'sample','cluster'};
    j = find(strncmpi(start,startNames,length(start)));
    if length(j) > 1
        error(message('litekmeans:AmbiguousStart', start));
    elseif isempty(j)
        error(message('litekmeans:UnknownStart', start));
    elseif isempty(k)
        error('litekmeans:MissingK', ...
            'You must specify the number of clusters, K.');
    end
    if j == 2
        if floor(.1*n) < 5*k
            j = 1;
        end
    end
    start = startNames{j};
elseif isnumeric(start)
    if size(start,2) == p
        center = start;
    elseif (size(start,2) == 1 || size(start,1) == 1)
        center = X(start,:);
    else
        error('litekmeans:MisshapedStart', ...
            'The ''Start'' matrix must have the same number of columns as X.');
    end
    if isempty(k)
        k = size(center,1);
    elseif (k ~= size(center,1))
        error('litekmeans:MisshapedStart', ...
            'The ''Start'' matrix must have K rows.');
    end
    start = 'numeric';
else
    error('litekmeans:InvalidStart', ...
        'The ''Start'' parameter value must be a string or a numeric matrix or array.');
end

% The maximum iteration number is default 100
if isempty(maxit)
    maxit = 100;
end

% The maximum iteration number for preliminary clustering phase on random
% 10% subsamples is default 10 
if isempty(clustermaxit)
    clustermaxit = 10;
end


% Assume one replicate
if isempty(reps) || ~isempty(center)
    reps = 1;
end

if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 0 && (round(k)==k))
    error('litekmeans:InvalidK', ...
        'X must be a positive integer value.');
elseif n < k
    error('litekmeans:TooManyClusters', ...
        'X must have more rows than the number of clusters.');
end


bestlabel = [];
sumD = zeros(1,k);
bCon = false;

for t=1:reps
    switch start
        case 'sample'
            center = X(randsample(n,k),:);
        case 'cluster'
            Xsubset = X(randsample(n,floor(.1*n)),:);
            [dump, center] = litekmeans(Xsubset, k, varargin{:}, 'start','sample', 'replicates',1 ,'MaxIter',clustermaxit);
        case 'numeric'
    end
    
    last = 0;label=1;
    it=0;
    
    switch distance
        case 'sqeuclidean'
            while any(label ~= last) && it<maxit
                last = label;
                
                bb = full(sum(center.*center,2)');
                ab = full(X*center');
                D = bb(ones(1,n),:) - 2*ab;
                
                [val,label] = min(D,[],2); % assign samples to the nearest centers
                ll = unique(label);
                if length(ll) < k
                    %disp([num2str(k-length(ll)),' clusters dropped at iter ',num2str(it)]);
                    missCluster = 1:k;
                    missCluster(ll) = [];
                    missNum = length(missCluster);
                    
                    aa = sum(X.*X,2);
                    val = aa + val;
                    [dump,idx] = sort(val,1,'descend');
                    label(idx(1:missNum)) = missCluster;
                end
                E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
                center = full((E*spdiags(1./sum(E,1)',0,k,k))'*X);    % compute center of each cluster
                it=it+1;
            end
            if it<maxit
                bCon = true;
            end
            if isempty(bestlabel)
                bestlabel = label;
                bestcenter = center;
                if reps>1
                    if it>=maxit
                        aa = full(sum(X.*X,2));
                        bb = full(sum(center.*center,2));
                        ab = full(X*center');
                        D = bsxfun(@plus,aa,bb') - 2*ab;
                        D(D<0) = 0;
                    else
                        aa = full(sum(X.*X,2));
                        D = aa(:,ones(1,k)) + D;
                        D(D<0) = 0;
                    end
                    D = sqrt(D);
                    for j = 1:k
                        sumD(j) = sum(D(label==j,j));
                    end
                    bestsumD = sumD;
                    bestD = D;
                end
            else
                if it>=maxit
                    aa = full(sum(X.*X,2));
                    bb = full(sum(center.*center,2));
                    ab = full(X*center');
                    D = bsxfun(@plus,aa,bb') - 2*ab;
                    D(D<0) = 0;
                else
                    aa = full(sum(X.*X,2));
                    D = aa(:,ones(1,k)) + D;
                    D(D<0) = 0;
                end
                D = sqrt(D);
                for j = 1:k
                    sumD(j) = sum(D(label==j,j));
                end
                if sum(sumD) < sum(bestsumD)
                    bestlabel = label;
                    bestcenter = center;
                    bestsumD = sumD;
                    bestD = D;
                end
            end
        case 'cosine'
            while any(label ~= last) && it<maxit
                last = label;
                W=full(X*center');
                [val,label] = max(W,[],2); % assign samples to the nearest centers
                ll = unique(label);
                if length(ll) < k
                    missCluster = 1:k;
                    missCluster(ll) = [];
                    missNum = length(missCluster);
                    [dump,idx] = sort(val);
                    label(idx(1:missNum)) = missCluster;
                end
                E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
                center = full((E*spdiags(1./sum(E,1)',0,k,k))'*X);    % compute center of each cluster
                centernorm = sqrt(sum(center.^2, 2));
                center = center ./ centernorm(:,ones(1,p));
                it=it+1;
            end
            if it<maxit
                bCon = true;
            end
            if isempty(bestlabel)
                bestlabel = label;
                bestcenter = center;
                if reps>1
                    if any(label ~= last)
                        W=full(X*center');
                    end
                    D = 1-W;
                    for j = 1:k
                        sumD(j) = sum(D(label==j,j));
                    end
                    bestsumD = sumD;
                    bestD = D;
                end
            else
                if any(label ~= last)
                    W=full(X*center');
                end
                D = 1-W;
                for j = 1:k
                    sumD(j) = sum(D(label==j,j));
                end
                if sum(sumD) < sum(bestsumD)
                    bestlabel = label;
                    bestcenter = center;
                    bestsumD = sumD;
                    bestD = D;
                end
            end
    end
end

label = bestlabel;
center = bestcenter;
if reps>1
    sumD = bestsumD;
    D = bestD;
elseif nargout > 3
    switch distance
        case 'sqeuclidean'
            if it>=maxit
                aa = full(sum(X.*X,2));
                bb = full(sum(center.*center,2));
                ab = full(X*center');
                D = bsxfun(@plus,aa,bb') - 2*ab;
                D(D<0) = 0;
            else
                aa = full(sum(X.*X,2));
                D = aa(:,ones(1,k)) + D;
                D(D<0) = 0;
            end
            D = sqrt(D);
        case 'cosine'
            if it>=maxit
                W=full(X*center');
            end
            D = 1-W;
    end
    for j = 1:k
        sumD(j) = sum(D(label==j,j));
    end
end
end



function [eid,emsg,varargout]=getargs(pnames,dflts,varargin)

% Initialize some variables
emsg = '';
eid = '';
nparams = length(pnames);
varargout = dflts;
unrecog = {};
nargs = length(varargin);

% Must have name/value pairs
if mod(nargs,2)~=0
    eid = 'WrongNumberArgs';
    emsg = 'Wrong number of arguments.';
else
    % Process name/value pairs
    for j=1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            eid = 'BadParamName';
            emsg = 'Parameter name must be text.';
            break;
        end
        i = strcmpi(pname,pnames);
        i = find(i);
        if isempty(i)
            % if they've asked to get back unrecognized names/values, add this
            % one to the list
            if nargout > nparams+2
                unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
                % otherwise, it's an error
            else
                eid = 'BadParamName';
                emsg = sprintf('Invalid parameter name:  %s.',pname);
                break;
            end
        elseif length(i)>1
            eid = 'BadParamName';
            emsg = sprintf('Ambiguous parameter name:  %s.',pname);
            break;
        else
            varargout{i} = varargin{j+1};
        end
    end
end

varargout{nparams+1} = unrecog;
end

% compute squared Euclidean distance
% ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
function d = L2_distance_1(a,b)
% a,b: two matrices. each column is a data
% d:   distance matrix of a and b



if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end

aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;

d = real(d);
d = max(d,0);

% % force 0 on the diagonal? 
% if (df==1)
%   d = d.*(1-eye(size(d)));
% end
end

function [x ft] = EProjSimplex_new(v, k)

%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=1
%

if nargin < 2
    k = 1;
end;

ft=1;
n = length(v); %O(n)

v0 = v-mean(v) + k/n; %O(n)
%vmax = max(v0);
vmin = min(v0); %O(n)
if vmin < 0
    f = 1; %O(1)
    lambda_m = 0; %O(1)
    while abs(f) > 10^-10 %O(n)
        v1 = v0 - lambda_m; %O(n)
        posidx = v1>0; %O(n)
        npos = sum(posidx); %O(n)
        g = -npos; %O(1)
        f = sum(v1(posidx)) - k; %O(n)
        lambda_m = lambda_m - f/g; %O(1)
        ft=ft+1; %O(1)
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0); %O(n)

else
    x = v0;
end;
%算法的时间复杂度是O(n)

end

function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)

if nargin < 2
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1)
    c = size(A,1);
end;

if nargin < 3
    isMax = 1;
    isSym = 1;
end;

if nargin < 4
    isSym = 1;
end;

if isSym == 1
    A = max(A,A');
end;
A(isnan(A)) = 0;
A(isinf(A)) = 0;
[v d] = eig(A); %v特征向量 d特征值
d = diag(d); %d本来是矩阵，取成向量
%d = real(d);
if isMax == 0
    [d1, idx] = sort(d);
else
    [d1, idx] = sort(d,'descend');
end;

idx1 = idx(1:c);
eigval = d(idx1); %前c个最大特征值
eigvec = v(:,idx1); %对应前c个最大特征值的特征向量

eigval_full = d(idx);

end

function [y1,y2,SS,U,V,cs] = coclustering_bipartite_fast1(A, c, NITER, islocal)

if nargin < 4
    islocal = 1;
end

if nargin < 3
    NITER = 30;
end

zr = 10e-11;
lambda = 0.1;

[n,m] = size(A);
onen = 1/n*ones(n,1);
onem = 1/m*ones(m,1);

A = sparse(A);
a1 = sum(A,2);
D1a = spdiags(1./sqrt(a1),0,n,n); 
a2 = sum(A,1);
D2a = spdiags(1./sqrt(a2'),0,m,m); 
A1 = D1a*A*D2a;

SS2 = A1'*A1; 
SS2 = full(SS2);


[V, ev0, ev]=eig1(SS2,m); %3-31 9：41  
aa = abs(ev); aa(aa>1-zr)=1-eps;
ad1 = aa(2:end)./aa(1:end-1);
ad1(ad1<0.15)=0; ad1 = ad1-eps*(1:m-1)'; ad1(1)=1;
ad1 = 1 - ad1;
[scores, cs] = sort(ad1,'descend');
cs = [cs, scores];

if nargin == 1
    c = cs(1);
end

V = V(:,1:c); 
U=(A1*V)./(ones(n,1)*sqrt(ev0(1:c)'));
U = sqrt(2)/2*U; V = sqrt(2)/2*V;  

a(:,1) = ev;
A = full(A); 



idxa = cell(n,1);
for i=1:n
    if islocal == 1
        idxa0 = find(A(i,:)>0);
    else
        idxa0 = 1:m;
    end
    idxa{i} = idxa0; 
end


idxam = cell(m,1);
for i=1:m
    if islocal == 1
        idxa0 = find(A(:,i)>0);
    else
        idxa0 = 1:n;
    end
    idxam{i} = idxa0; 
end

D1 = 1; D2 = 10;
for iter = 1:NITER
    
    U1 = D1*U;
    V1 = D2*V;
    dist = L2_distance_1(U1',V1');  % only local distances need to be computed. speed will be increased using C

    S = zeros(n,m);
    for i=1:n
        idxa0 = idxa{i};
        ai = A(i,idxa0);
        di = dist(i,idxa0);
        ad = (ai-0.5*lambda*di); 
        %S(i,idxa0) = EProjSimplex_new(ad);
        
        nn = length(ad);
        %v0 = ad-mean(ad) + 1/nn;
        v0 = ad-sum(ad)/nn + 1/nn;
        vmin = min(v0);
        if vmin < 0
            lambda_m = 0;
            while 1
                v1 = v0 - lambda_m;
                %posidx = v1>0; npos = sum(posidx);
                posidx = find(v1>0); npos = length(posidx);
                g = -npos;
                f = sum(v1(posidx)) - 1;
                if abs(f) < 10^-6
                    break;
                end
                lambda_m = lambda_m - f/g;
            end
            vv = max(v1,0);
            S(i,idxa0) = vv;
        else
            S(i,idxa0) = v0;
        end
    end
    
    
    %Sm = sparse(m,n);
    
    Sm = zeros(m,n);
    for i=1:m
        idxa0 = idxam{i};
        ai = A(idxa0,i);
        di = dist(idxa0,i);
        ad = (ai-0.5*lambda*di);
        Sm(i,idxa0) = EProjSimplex_new(ad);
    end
    
    S = sparse(S);
    Sm = sparse(Sm);    
    SS = (S+Sm')/2;
    %SS = sparse(SS);
    d1 = sum(SS,2);
    D1 = spdiags(1./sqrt(d1),0,n,n);
    d2 = sum(SS,1);
    D2 = spdiags(1./sqrt(d2'),0,m,m);
    SS1 = D1*SS*D2;
    
    SS2 = SS1'*SS1;
    SS2 = full(SS2);
    [V, ev0, ev]=eig1(SS2,c);
    U=(SS1*V)./(ones(n,1)*sqrt(ev0'));
    U = sqrt(2)/2*U; V = sqrt(2)/2*V;
    

    
    a(:,iter+1) = ev;
    U_old = U;
    V_old = V;
    
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c)); 
    if fn1 < c-0.0000001
        lambda = 2*lambda;
    elseif fn2 > c-0.0000001
        lambda = lambda/2;   U = U_old; V = V_old; 
    else
        break;
    end
end

SS0=sparse(n+m,n+m); SS0(1:n,n+1:end)=SS; SS0(n+1:end,1:n)=SS';  
[clusternum, y]=graphconncomp(SS0);

y1=y(1:n)'; 
y2=y(n+1:end)'; 


end

