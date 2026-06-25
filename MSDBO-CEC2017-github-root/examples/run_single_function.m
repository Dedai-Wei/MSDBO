clear;
clc;

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'benchmarks'));

funcNum = 1;
popSize = 30;
maxIter = 300;
dim = 30;
lb = -100;
ub = 100;

rng(2026, 'twister');
fobj = @(x) cec2017_eval(x, funcNum);
[bestScore, bestPosition, curve] = MSDBO(popSize, maxIter, lb, ub, dim, fobj);

fprintf('Function: F%02d\n', funcNum);
fprintf('Best score: %.6e\n', bestScore);
fprintf('Best position norm: %.6e\n', norm(bestPosition));

figure;
semilogy(curve, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Best fitness');
title(sprintf('MSDBO on F%02d', funcNum));
grid on;
