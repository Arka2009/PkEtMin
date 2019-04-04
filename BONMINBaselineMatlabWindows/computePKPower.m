% Compute the power of a single phase
function y = computePKPower(c,benchid,AP,BP,LLIM,M,N)
    % if ((ci > N) || ci < LLIM(benchid))
    %     error('computePower:PKP','Allocaction of %d core to bench-%d is illegal',ci,benchid)
    % end
    z = zeros(M,1);
    for ph=1:M
        a     = AP(benchid(ph));
        b     = BP(benchid(ph));
        z(ph) = a*c(ph) + b;
    end
    y = max(z);
end