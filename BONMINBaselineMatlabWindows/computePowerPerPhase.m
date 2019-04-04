% Compute the power of a single phase
function y = computePowerPerPhase(ci,benchid,AP,BP,LLIM,M,N)
    % if ((ci > N) || ci < LLIM(benchid))
    %     error('computePower:PKP','Allocaction of %d core to bench-%d is illegal',ci,benchid)
    % end
    a = AP(benchid);
    b = BP(benchid);
    y = a*ci + b;
end