% Stress test integer-based matrix multiply of Ootomo et al. [1] and
% demonstrate scenarios in which it may be highly inaccurate by using
% a minimal example of an inner product of two 2-element vectors.
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
% [2] Error Analysis of Floating-Point Matrix Multiplication Computed
%     via Low-Precision Integer Arithmetic. Ahmad Abdelfattah,
%     Jack Dongarra, Massimiliano Fasi, Mantas Mikaitis, and
%     Francoise Tisseur. arXiv:2505.XXXXX [math.NA]. June, 2025.

rng("default")

clear t_vals s_vals
global t_vals s_vals
t_vals = 0:100;  % Parameter that controls the spread of ranges.
s_vals = 1:30;   % Number of splits.

fwd_err_gemmi_ooy = zeros(length(t_vals), length(s_vals));
fwd_err_gemmi_uoi = zeros(length(t_vals), length(s_vals));
fwd_err_binary64 = zeros(length(t_vals), 1);

% ooy24 uses bit splitting and floating-point accumulation.
algin_ooy.split = 'b';
algin_ooy.acc = 'f';

% uoi24 uses round-to-nearest and integer accumulation.
algin_uoi.split = 'n';
algin_uoi.acc = 'i';

for i = 1:length(t_vals)
  t = t_vals(i);
  a = randn() * 2 .^ [-t, 0];
  b = randn() * 2 .^ [t; 0];
  ref = sym(a) * sym(b);  % Reference solution.
  fwd_err_binary64(i) = abs(a * b - ref) / abs(ref); % Binary64 solution.
  for j = 1:length(s_vals)
    s = s_vals(j);
    fwd_err_gemmi_ooy(i,j) =...
        abs(gemmi(a, b, s, s, algin_ooy) - ref) / abs(ref);
    fwd_err_gemmi_uoi(i,j) =...
        abs(gemmi(a, b, s, s, algin_uoi) - ref) / abs(ref);
  end
end

% Find minimum and maximum values for color map.
max_exponent = log10(max(abs([fwd_err_gemmi_ooy(:)])));
min_exponent = log10(min(abs([fwd_err_gemmi_ooy(:)])));

print_results(fwd_err_binary64, './data/fwd_err_binary64.dat')
print_results(fwd_err_gemmi_ooy, './data/fwd_err_gemmi_ooy.dat')
print_results(fwd_err_gemmi_uoi, './data/fwd_err_gemmi_uoi.dat')

% Plot results.
subplot(1, 2, 1)
imagesc(log10(fwd_err_gemmi_ooy))
title('Ootomo et al.')
xlabel("s");
ylabel("t");
shading interp;
caxis manual
caxis([min_exponent max_exponent]);

subplot(1, 2, 2)
imagesc(log10(fwd_err_gemmi_uoi))
title('Uchino et al.')
xlabel("s");
ylabel("t");
shading interp;
caxis manual
caxis([min_exponent max_exponent]);

colorbar

% Export result.
function print_results(variable, filename)
  global t_vals s_vals
  output_file = fopen(filename, 'w');
  fprintf(output_file, 't s magnitude\n');
  for j = 1:size(variable, 2)
    for i = 1:size(variable, 1)
      fprintf(output_file, '%3i %3i %.20f\n',...
          t_vals(i), s_vals(j), log10(variable(i, j)));
    end
  end
end