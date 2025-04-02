
function [preY,objres,sum_Z] = HMIC(X,c,m,alpha,beta,rho,lp_param,seed,kernel_type,norm_type,opts)
rand('seed',seed)
if (~exist('opts','var'))
   opts = [];
end

Distance = 'sqEuclidean';  %(the default)
if isfield(opts,'Distance')
    Distance = opts.Distance;
end

%For quadratic programming (QP) problem, there are two options to obtain
%the solution (i.e. 'SimplexQP_acc' or 'quadprog').
QP_options = 'quadprog';
%QP_options = 'SimplexQP_acc';


View = length(X); %The number of views
[n,~] = size(X{1}); %The number of samples
weight_vector = ones(1,View)';
%Maximum and minimum normalization
XX = [];
for i = 1:View

    X{i} = ( X{i}-repmat(min(X{i}),n,1) ) ./repmat( max(X{i})-min(X{i}),n,1);
    X{i}( isnan(X{i}) ) = 1;
    XX = [XX X{i}];
    d(i) = size(X{i},2); %The number of features in the i-th view
end


    subfea = XX;


% rand('twister',seed);
[~, AA] = litekmeans(subfea,m,'MaxIter', 20,'Replicates',2,'Distance',Distance); %Calculate the distance based on the input 'Distance'

%Split
temp = 0;
for ia = 1:View
    A{ia} = AA(:,1 + temp : d(ia) + temp);
    temp = temp + d(ia);
end

%Initialization and Optimization 
Itermax = 20;
IterMaxP = 50;

delta = 1 / View * ones(View,1);

for v = 1:View
    G{v} = zeros(n,m);
    W{v} = zeros(n,m);
end
w = zeros(n*m*View,1); g = zeros(n*m*View,1);
sX = [n,m,View];
%% 为每个视图创建options，其中包括每个视图的核函数类型等，用于后续的核生成。
for i=1:View
    options(i).KernelType = kernel_type;
    options(i).t = optSigma(X{i});
    options(i).d = 4;
end
%KNN  init U
opts1.disp = 0;
U = cell(1,View);
K = cell(1,View);
L = cell(1,View);
for iv = 1:View
    epsilon = 1e-6;  % 小的正则化项

    K{iv} = constructKernel(X{iv},X{iv},options(iv));
    D = diag(sum(K{iv},1));
    L{iv} = sqrt(inv(D))*K{iv}*sqrt(inv(D));  
    L{iv}=(L{iv}+L{iv})/2;

    
    [U{iv} E] = eigs(L{iv},m,'LA',opts1);  
end


%% Init P
P = delta(1) * U{1}; 
for iv = 2:View
    P = P + delta(iv) * U{iv};
end

res = zeros(Itermax + 1,1); %Obj_value
res(1) = obj_value(X,A,U,P,View,alpha,beta,delta);
deltaArray = zeros(View,Itermax + 1);
deltaArray(:,1) = delta;
fprintf('iter = 0, obj_value = %f\n',res(1))
for i = 1:Itermax
    %Update P    
    disp([num2str(i),'-th iteration...']);
    tic1 = tic;
    sum_U = sparse(n,m);
    for i1 = 1:View
        sum_U = sum_U + delta(i1) * U{i1};
    end
    sum_U(find(isnan(sum_U)==1)) = 0;

    [~,~,P,~,~,~] = coclustering_bipartite_fast1(sum_U, c, IterMaxP);
    

    %% 大改
    %% Update Uv
    U = updateU(X,A,alpha,beta,delta,P,U,W,G,rho);
    
    %% Update G
    U_tensor = cat(3, U{:,:});
    W_tensor = cat(3, W{:,:});

    u = U_tensor(:);
    w = W_tensor(:);
    % 给U加一个nan 和inf = 0 和 1 的代码
    
    if(norm_type == 'SP')
        [g, objV] = wshrinkObj_weight_lp(1/rho*w - u,weight_vector*rho/2,sX,0,3,lp_param);
    elseif(norm_type == 'Nuclear')
        [g, objV] = wshrinkObj_weight(1/rho*w - u,weight_vector*rho/2,sX,0,3);
    end
    G_tensor = reshape(g, sX);
    for k=1:View
        G{k} = G_tensor(:,:,k);
    end
    %% Update w
    w = w + rho*(u - g);
    W_tensor = reshape(w, sX);
    for k=1:View
        W{k} =  W_tensor(:,:,View);
    end
    
    %Update delta
    for iv = 1:View
        U{iv} = sparse(U{iv});
    end
    UU = sparse(n * m,View);
    for i2 = 1:View
        UU(:,i2) = reshape(U{i2},[n*m 1]);
    end
    newU = UU'*UU;
    p = reshape(P,[n*m 1]);
    s = 2*UU'*p;
    %QP_options
    switch lower(QP_options)
        case {lower('SimplexQP_acc')}
            delta = SimplexQP_acc(newU, s);
            
        case {lower('quadprog')}
            options = optimset( 'Algorithm','interior-point-convex','Display','off');
            delta = quadprog(2*newU,-s,[],[],ones(1,View),1,zeros(View,1),ones(View,1),[],options);            
    end 
    deltaArray(:,i + 1) = delta;
    
    
    
    
    %Calculate the objective function value
    res(i+1) = obj_value(X,A,U,P,View,alpha,beta,delta);
    if abs(res(i+1) - res(i)) < 1e-5 || norm(deltaArray(:,i + 1) - deltaArray(:,i),2) < 1e-5
        break
    end   
    objres = res;
    toc(tic1);
end


sum_Z = 0;
for i1 = 1:View
    sum_Z = sum_Z + delta(i1) * U{i1};
end



preY = coclustering_bipartite_fast1(sum_Z, c, IterMaxP);



end

