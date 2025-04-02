function norm_21 = l21_norm(U)
    norm_21 = sum(sqrt(sum(U.^2, 2))); % 对每一行计算2范数，并对结果求和
end