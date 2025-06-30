%DISTANCETOLASTBIT1 Number of bits effectively used in a significand.
%
% Compute the distance between the most significant bit set to 1,
% inclusive, and the least significant bit set to 1, inclusive,
% in the significand. For the normalized floating point values the implicit
% bit is included; assumes the inputs are encoded in IEEE 754 binary64
% format encoding.
%
% References:
%
% [1] Analysis of Floating-Point Matrix Multiplication Computed
%     via Integer Arithmetic. Ahmad Abdelfattah, Jack Dongarra,
%     Massimiliano Fasi, Mantas Mikaitis, and Francoise Tisseur.
%     arXiv:2506.11277 [math.NA]. June, 2025.

function A = distanceToLastBit1(A)
    [n, k] = size(A);
    for i = 1 : n
        for j = 1 : k
            if A(i, j) ~= 0
                [~, exponent] = log2(A(i, j));
                exponent = exponent - 1;
                bit_pattern = string(dec2bin(typecast(A(i, j), 'uint64')));
                orig_length = strlength(bit_pattern);
                temp = strip(bit_pattern, "0");
                trail_zeros = min(52, orig_length - strlength(temp));
                A(i, j) = 53 - trail_zeros;
                if (exponent < -1022)
                    A(i, j) = A(i, j) - (-1022 - exponent);
                end
            end
        end
    end
end
