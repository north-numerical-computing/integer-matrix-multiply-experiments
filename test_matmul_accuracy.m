% Reproduce an experimental testbench of Ootomo et al. [1] for measuring
% the accuracy of integer-based matrix multiplication.
%
% Requirements:
%
%   gemmi: https://github.com/north-numerical-computing/gemmi
%
% References:
%
% [1] H. Ootomo, K. Ozaki, and R. Yokota. DGEMM on integer matrix
%     multiplication unit. Int. J. High. Perf. Comput. Appl. 2024.
%
% [2] Analysis of Floating-Point Matrix Multiplication Computed
%     via Integer Arithmetic. Ahmad Abdelfattah, Jack Dongarra,
%     Massimiliano Fasi, Mantas Mikaitis, and Francoise Tisseur.
%     arXiv:2506.11277 [math.NA]. June, 2025.

clear all;
rng(0)

% Set up gemmi to use splitting based on bit-level operations.
alg.split = 'b';

% Different split settings for gemmi.
splits = {2, 4, 6, 8, 10};

% Matrix dimensions: A is m x n, B is n x q
m = 10;
q = 10;
nlist = floor(logspace(1, 5, 20));

% Parameter from experiments of [1] that controls the wideness of exponent
% distribution of rows of A and columns of B.
phis = {8, 13};

tiledlayout(1, 2);

for phi = phis
    for n = 1:length(nlist)
        % Generate the test matrices.
        A = (rand(m, nlist(n)) - 0.5) .* exp(phi{1} * randn(m, nlist(n)));
        B = (rand(nlist(n), q) - 0.5) .* exp(phi{1} * randn(nlist(n), q));

        % Compute a reference result with 32 digits using vpa.
        Ctrue = double(vpa(A) * vpa(B));

        % Compute binary64 errors.
        Cdouble = A * B;
        err_double(n) = max(...
            abs(Cdouble - Ctrue) ./ abs(Ctrue), [], 'all');

        for s = 1:length(splits)
            C = gemmi(A, B, splits{s}, splits{s}, alg);

            % Maximum forward error
            err_oz(s, n) = max(abs(C - Ctrue) ./ abs(Ctrue), [], 'all');
        end
    end

    nexttile
    loglog(nlist, err_oz(1, :), '-o');
    hold on
    loglog(nlist, err_oz(2, :), '--o');
    loglog(nlist, err_oz(3, :), '-*');
    loglog(nlist, err_oz(4, :), '--*');
    loglog(nlist, err_oz(5, :), '-x');
    loglog(nlist, err_double, '--x');
    legend('Split:2', 'Split:4', 'Split:6',...
        'Split:8','Split 10','Standard binary64');
    hold off

    % Output various results to .dat files.
    for i = 1:length(phi)
        filename = strcat('data/test_matmul_berr_accuracy_phi',...
            num2str(phi{i}), '.dat');
        fileID = fopen(filename, 'w');
        fprintf(fileID,...
           'n split2 split4 split6 split8 split10 standard-binary64 \n');
        for j = 1:length(nlist)
            fprintf(fileID,'%d %e %e %e %e %e %e \n',...
                nlist(j), err_oz(1, j), err_oz(2, j), err_oz(3, j),...
                err_oz(4, j), err_oz(5, j), err_double(j));
        end
    end
end
