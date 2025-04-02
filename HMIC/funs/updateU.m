function U = updateU(X,A,alpha,beta,delta,P,initU,W,G,rho)
%For quadratic programming (QP) problem, there are two options to obtain
%the solution (i.e. 'SimplexQP_acc' or 'quadprog').
% QP_options = 'quadprog';
%QP_options = 'SimplexQP_acc';

[m,~] = size(A{1}); %The number of anchors
n = size(X{1},1); %The number of samples
View = length(X); %The number of views
if nargin < 4 %Initialize Z
    beta = 0;
    delta = zeros(View,1);
    P = zeros(n,m);
    for j = 1:View
        initU{j} = P;
    end
end
% options = optimset( 'Algorithm','interior-point-convex',' Display','off');
U = cell(1,View);
for i = 1:View
    U{i} = zeros(n, m); % 或者根据你的需求初始化矩阵
end
%協同正則化部分

for i = 1:View
%     sum_UU = 0;
    sum_Ui = 0;
    sum_UU = 0;
    for j = 1:View
        if i == j
            continue;
        end
        sum_UU = sum_UU + initU{j}'*initU{j};
        sum_Ui =sum_Ui + delta(j)*U{j};
    end

    
    
    %% 偏导推导
%     U1 = A{i}*A{i}' - alpha*sum_UU +( beta * delta(i)^2 + rho) * eye(m);
%     U1 = inv(U1);
% %     U2 = X{i}*A{i}' - beta*delta(i)*(sum_Ui - P) + W{i} - rho*G{i};
%     U2 = X{i}*A{i}' - beta*delta(i)*(sum_Ui - P);
       
    U1 = 2*A{i}*A{i}' + (2*beta*delta(i)^2 + rho)* eye(m) - 2*alpha * sum_UU ;
    U1 = inv(U1);
    U2 = 2*X{i}*A{i}' - 2*beta*delta(i)*(sum_Ui - P) - W{i} +rho*G{i} ;
    
   %% 源代码QP
%     H = L{i} + alpha*sum_UU + (beta * delta(i)^2)*eye(n) +A{i}*A{i}';
%     H = (H+H')/2;
%     B = A{i};
% 
%     delU = zeros(n,m);
%     for i1 = 1:View
%         if i1 ~= i
%             delU = delU + delta(i1) * initU{i1};
%         end
%     end
%     Uv = zeros(n,m);
%     for ji = 1:m
% %         size(B(ji,:))
% %         size(X{i}')
% %         size(delU(:,ji)')
% %         
% %         tmp1 = 2*B(ji,:)*X{i}';
% %         tmp2 = 2 * beta * delta(i) * (delU(:,ji)' - P(:,ji)');
%         ff = -2*B(ji,:)*X{i}'  + 2 * beta * delta(i) * (delU(:,ji)' - P(:,ji)');
%         %QP_options
%         switch lower(QP_options)
%             case {lower('SimplexQP_acc')}
%                 Uv(:,ji) = SimplexQP_acc(H / 2,-ff');                
%             case {lower('quadprog')}
%                 Uv(:,ji) = quadprog(H,ff',[],[],ones(1,n),1,zeros(n,1),ones(n,1),[],options);
% 
%         end        
%     end
    U{i}=U2*U1;


end
