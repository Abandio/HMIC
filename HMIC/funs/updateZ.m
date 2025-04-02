function Z = updateZ(X,A,alpha,beta,delta,P,initZ)
%For quadratic programming (QP) problem, there are two options to obtain
%the solution (i.e. 'SimplexQP_acc' or 'quadprog').
QP_options = 'quadprog';
%QP_options = 'SimplexQP_acc';

[m,~] = size(A{1}); %The number of anchors
n = size(X{1},1); %The number of samples
View = length(X); %The number of views
if nargin < 4 %Initialize Z
    beta = 0;
    delta = zeros(View,1);
    P = zeros(n,m);
    for j = 1:View
        initZ{j} = P;
    end
end
options = optimset( 'Algorithm','interior-point-convex','Display','off');
Z = cell(1,View);
for i = 1:View
%     H = 2 * (alpha + beta * delta(i)^2) * eye(m) + 2*A{i}*A{i}';
    H = 2 * (A{i}*A{i}' +(beta * delta(i)^2)*eye(m) - 2*( P'*P))
    H = (H+H')/2;
    B = X{i}';
    delZ = zeros(n,m);
    for i1 = 1:View
        if i1 ~= i
            delZ = delZ + delta(i1) * initZ{i1};
        end
    end
    Zv = zeros(size(P,2),n);
    parfor ji = 1:n
        ff = -2*B(:,ji)'*A{i}' + 2 * beta * delta(i) * (delZ(ji,:) - P(ji,:));
        %QP_options
        switch lower(QP_options)
            case {lower('SimplexQP_acc')}
                Zv(:,ji) = SimplexQP_acc(H / 2,-ff');                
            case {lower('quadprog')}
                Zv(:,ji) = quadprog(H,ff',[],[],ones(1,m),1,zeros(m,1),ones(m,1),[],options);
        end        
    end
    Z{i}=Zv';
end

end

