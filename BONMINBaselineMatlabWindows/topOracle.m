% Load the benchmark characteristics 
function xyz = topOracle(M,suffix)
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
        D       = w(1);
        benchid = w(2:M+1) + 1; % MATLAB follows 1-indexed convention
        
        % disp(benchid);
        maxtime = 1.05*(10^(M/10));
        tic;
        optx    = computeOptimalOracle(D,benchid,AET,BET,AP,BP,LLIM,M,N);
        elapsed = toc;

        
        % Dump the output
        for i=1:M
            fprintf(fild3,'%d,',optx(i));
        end
        fprintf(fild3,'%f,',optx(M+1));
        if computeExecTime((optx(1:M))',benchid,AET,BET,LLIM,M,N) <= D
            fprintf(fild3,'passed,');
        else
            fprintf(fild3,'failed,');
        end
        fprintf(fild3,'%f\n',elapsed);
        disp(optx);

        % Break if you have found a point with a high execution time
        if elapsed >= maxtime
            fprintf('Found a large execution time, exiting')
            break
        end
    end
    fprintf('Completed Experiments with M = %d phases, N = %d cores, elapsedTime = %f\n',M,N,elapsed);
    fclose(fild3);
end
