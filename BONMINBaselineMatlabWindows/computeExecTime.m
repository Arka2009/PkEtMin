% Execution Time
function y = computeExecTime(c,benchid,AET,BET,LLIM,M,N)
    % y = a*x + b
    y = 0;
    for i=1:M
        a = AET(benchid(i));
        b = BET(benchid(i));
        % if ((c(i) > N) || c(i) < LLIM(benchid(i)))
        %     error('computeExecTime:ET','Allocaction to to %d-phase is illegal',i)
        % end  
        y = y + 1/((a*log(c(i))+b));
    end
end