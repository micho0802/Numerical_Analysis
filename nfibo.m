
function [x] = nfibo(n)
    for i = 1:n
        if i == 1;
            x(i) = 0;
        elseif i == 2;
            x(i) = 1;
        else
            x(i) = x(i-1) + x(i-2);
        end
    end
