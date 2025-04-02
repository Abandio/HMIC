%  min  x'*A*x - x'*b
%  s.t. x'1=1, x>=0

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