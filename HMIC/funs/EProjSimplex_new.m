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
