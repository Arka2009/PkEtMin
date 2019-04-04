% Load the benchmark characteristics 
function xyz = topOracleDual(M,suffix)
    % Number of cores is fixed to 16
    N = 16;

    load('benchParams.mat');
    % if (exist AET  ~= 1 || ...
    %     exist BET  ~= 1 || ...
    %     exist AP   ~= 1 || ...
    %     exist BP   ~= 1 || ...
    %     exist LLIM ~= 1)
    %     error('All variables could not be loaded');
    % end

    % Read the workload from CSV file
    fil2 = sprintf('workloads-%s/wkld_%d.csv',suffix,M);
    M2   = csvread(fil2);

    % Dump Matlab Output
    fil3  = sprintf('workloads-%s/wkld_%d_matlab.out.csv',suffix,M);
    fild3 = fopen(fil3,'w');

    % Iterate through all the workloads
    for w2 = M2.'
        w       = w2';
        P       = w(1);
        benchid = w(2:M+1) + 1; % MATLAB follows 1-indexed convention
        
        % disp(benchid);
        tic;
        optx    = computeOptimalOracleDual(P,benchid,AET,BET,AP,BP,LLIM,M,N);
        elapsed = toc;

        
        % Dump the output
        for i=1:M
            fprintf(fild3,'%d,',optx(i));
        end
        fprintf(fild3,'%f,',computeExecTime(optx,benchid,AET,BET,LLIM,M,N));
        if computePKPower(optx',benchid,AET,BET,LLIM,M,N) <= P
            fprintf(fild3,'passed,');
        else
            fprintf(fild3,'failed,');
        end
        fprintf(fild3,'%f\n',elapsed);
        disp(optx);
    end
    fprintf('Completed Experiments with M = %d phases, N = %d cores, elapsedTime = %f\n',M,N,elapsed);
    fclose(fild3);
end