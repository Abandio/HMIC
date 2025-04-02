function res = obj_value(X,A,Z,P,View,alpha,beta,delta)
res = 0;
sum_Z = 0;
for i = 1:View
    res = res + norm(X{i}' - A{i}' * Z{i}','fro')^2 + alpha * norm(Z{i},'fro')^2;
    sum_Z = sum_Z + delta(i) * Z{i};
end
res = res + beta * norm(sum_Z - P,'fro')^2;
end