% Factorize A, solve a linear system, and output the accuracy, for
% an anymatrix matrix. Inside LU factorisation, integer-based matrix
% multiplication is utilised instead of standard binary64 arithmetic.
%
% Requirements:
%
%   gemmi: https://github.com/north-numerical-computing/gemmi
%   anymatrix: https://github.com/north-numerical-computing/anymatrix
%
% References:
%
% [1] Analysis of Floating-Point Matrix Multiplication Computed
%     via Integer Arithmetic. Ahmad Abdelfattah, Jack Dongarra,
%     Massimiliano Fasi, Mantas Mikaitis, and Francoise Tisseur.
%     arXiv:2506.11277 [math.NA]. June, 2025.

clear all;
rng(0);

% 1 - use integer-based GEMM, 0 - use standard binary64.
set_gemm_mode(1);

% Inputs to matrices that generate with a single input argument.
% In most cases it is one of the matrix dimensions controlling the size of
% the generated matrix. Matrices not following this pattern are excluded.
matrix_size = 500;

% Block size for the block LU factorisation.
block_size = 10;

rel_err = 0;
% Get IDs of all square dense matrices of real values. Triangular
% matrices are also excluded because they only require substitution
% when solving Ax=b.
matrix_ID = 'core/blockhouse';

% Try the generator with one scalar input arg.
try
A = anymatrix(matrix_ID, matrix_size);
catch
    return
end

b = rand(length(A), 1);

lastwarn('');
% Compute the reference solution in binary64.
x = A\b;
[warnMsg, warnId] = lastwarn;

% If matrix was generated as expected, and MATLAB could solve Ax=b
% without warnings about the large condition number, use this matrix.
if isempty(warnMsg) && ~isempty(A) && length(A) == matrix_size ...
        && ~ismember(inf, A) && ~ismember(-inf, A)...
        && ~ismember(nan, A) && ~ismember(inf, x) &&...
        ~ismember(-inf, x) && ~ismember(nan, x)
    fprintf("Solving Ax=b for matrix %s\n", matrix_ID);

    for split_A = 1:8
        for split_B = 1:8
            [x_comp, msplitsA, msplitsB] =...
                lu_block_pivoting_solve(...
                A, b, block_size, split_A, split_B);
            splits_needed_A(split_A, split_B) = max(msplitsA);
            splits_needed_B(split_A, split_B) = max(msplitsB);
    
            rel_err(split_A, split_B) =...
                norm(A * x_comp-b, inf)/...
                (eps * (norm(A, inf)*norm(x_comp, inf)+norm(b, inf)) *...
                matrix_size); % HPL error measurement.
        end
    end
end

% Output various results to .dat files.
filename = strcat('data/gaussian_IMMA_test_exhaustive_splits.dat');
fileID = fopen(filename, 'w');
fprintf(fileID, 'splitA splitB rel_error \n');
for split_A = 1:8
    for split_B = 1:8
    fprintf(fileID, '%d %d %e \n',...
        split_A, split_B, log(rel_err(split_A, split_B)));
    end
end


function OZAKI = set_gemm_mode(mode)
    global OZAKI;
    OZAKI = mode;
end

% Solve Ax=b by factorizing A with block Gaussian, with partial pivoting.
function [x, msplitsA, msplitsB] =...
    lu_block_pivoting_solve(A, b, r, splits_A, splits_B)

    temp_A = A;
    [L, U, piv, msplitsA, msplitsB] =...
        lu_block_pivoting(A, r, splits_A, splits_B);

    for k = 1:(length(A) - 1)
        temp = b(k);
        b(k) = b(piv(k));
        b(piv(k)) = temp;
    end

    x = U\(L\b);
end

% Block Gaussian LU factorization (Sec. 3.6.1 of MC4) with part. pivoting.
function [L, U, piv, msplitsA, msplitsB] =...
    lu_block_pivoting(A, r, splits_A, splits_B)

    global OZAKI;

    alg.split = 'b';

    n = length(A);
    N = n/r;

    L = zeros(n,n);
    U = zeros(n,n);

    gemms = 0;

    piv = [];

    for k=1:N
        rows = ((k - 1) * r + 1):n;
        columns = ((k - 1) * r + 1):((k - 1) * r + r);
        Ltemp = [];
        Utemp = [];
        [Ltemp, Utemp, pivpart] =...
            lu_rectangular_part_pivot(A(rows, columns));

        for i = ((k - 1) * r + 1):((k - 1) * r + length(pivpart))
            pivpart(i - (k - 1) * r) =...
                pivpart(i - (k - 1) * r) + (k - 1) * r;
            temp = A(i, :);
            A(i, :) = A(pivpart(i - (k - 1) * r), :);
            A(pivpart(i - (k - 1) * r), :) = temp;

            if (k >= 2)
                temp = L(i, :);
                L(i, :) = L(pivpart(i - (k - 1) * r), :);
                L(pivpart(i - (k - 1) * r), :) = temp;
            end
        end

        L(rows, columns) = Ltemp;
        U(columns, columns) = Utemp;

        piv = [piv, pivpart];

        for l = (k + 1):N
            U(columns, ((l - 1) * r + 1):((l - 1) * r + r)) =...
                L(columns, columns)\...
                A(columns, ((l - 1) * r + 1):((l - 1) * r + r));
        end

        if (k ~= N)
            if (OZAKI)
                A((k * r + 1):n, (k * r + 1):n) =...
                    A((k * r + 1):n, (k * r + 1):n) - ...
                    gemmi(L((k * r + 1):n, ((k - 1) * r + 1):(k * r)),...
                          U(((k - 1) * r + 1):(k * r), (k * r + 1):n), ...
                    splits_A, splits_B, alg);
                [~, tempA] =...
                    log2(L((k * r + 1):n, ((k - 1) * r + 1):(k * r)));
                [~, tempB] =...
                    log2(U(((k - 1) * r + 1):(k * r), (k * r + 1):n));
                tempA = tempA + 1;
                tempB = tempB + 1;
                gemms = gemms + 1;
                msplitsA(gemms) =...
                    ceil(max(max(max(tempA, [], 2)-...
                        tempA+distanceToLastBit1(...
                        L((k * r + 1):n, ((k - 1) * r + 1):(k * r))),...
                        [], 2))/7);
                msplitsB(gemms) =...
                    ceil(max(max(max(tempB, [], 1)-...
                    tempB + distanceToLastBit1(...
                    U(((k - 1) * r + 1):(k * r), (k * r + 1):n)),...
                    [], 1))/7);
            else
                A((k * r + 1):n, (k * r + 1):n) =...
                    A((k * r + 1):n, (k * r + 1):n) - ...
                    L((k * r + 1):n, ((k - 1) * r + 1):(k * r)) *...
                    U(((k - 1) * r + 1):(k * r), (k * r + 1):n);
            end
        end
    end
end

% Gaussian LU factorization with partial pivoting (Alg. 3.4.1 of MC4).
function [L, U, piv] = lu_simple_part_pivot(A)
    n = length(A);
    piv = [];

    for k = 1: (n-1)
        [~, mu] = max(abs(A(k:n, k)));
        mu = mu + (k - 1);
        piv(k) = mu;
        temp = A(k, :);
        A(k,:) = A(mu, :);
        A(mu, :) = temp;
        
        if (A(k, k) ~= 0)
            rows = (k + 1):n;
            A(rows, k) = A(rows, k) / A(k, k);
            A(rows, rows) = A(rows, rows) - A(rows, k) * A(k, rows);
        end
    end
    
    L = tril(A, -1) + eye(n, n);
    U = triu(A);
end

% Rectangular (tall and narrow) case of above.
function [L, U, piv] = lu_rectangular_part_pivot(A)
   [n,r] = size(A);
    if n == r
        [L, U, piv] = lu_simple_part_pivot(A);
        return;
    elseif n > r
        for k = 1:r
            [~, mu] = max(abs(A(k:n, k)));
            mu = mu + (k - 1);
            piv(k) = mu;
            temp = A(k, :);
            A(k, :) = A(mu, :);
            A(mu, :) = temp;

            if (A(k, k) ~= 0)
                rows = (k + 1):n;
                A(rows, k) = A(rows, k) / A(k, k);
                if (k<r)
                    columns = (k + 1):r;
                    A(rows, columns) = A(rows, columns) -...
                        A(rows, k) * A(k, columns);
                end
            end
        end
    end

    L = tril(A, -1) + eye(n, r);
    U = triu(A);
    U = U(1:r, 1:r);
end
