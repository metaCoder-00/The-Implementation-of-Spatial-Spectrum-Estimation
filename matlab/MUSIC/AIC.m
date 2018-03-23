function n = AIC(L, M, D)
%---n: the number of sources using AIC method---%
%---L: the number of snapshots------------------%
%---M: the number of array elements-------------% 
%---D: eigen values with descendant sort--------%
    t = inf;                
    for n = 1: M - 1
        num = 0;
        den = 1;
        for iter = n + 1: M
            num = num + D(iter);
            den = den * D(iter);
        end
        Lambda = (num/(M - n)) / den^(1/(M - n));
        etrpy = 2*L*(M - n)*log(Lambda) + 2*n*(2*M - n);    % AIC entropy
%-------find minimum loss entropy----------------%
        if etrpy > t
            break;
        else
            t = etrpy;
        end
    end
    n = n - 1;          % fix n
end