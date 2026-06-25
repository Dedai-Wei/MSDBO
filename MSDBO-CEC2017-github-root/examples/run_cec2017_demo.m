clear;
clc;

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'src'));
addpath(fullfile(rootDir, 'benchmarks'));

popSize = 30;
maxIter = 500;
dim = 30;
lb = -100;
ub = 100;
funcList = 1:30;

results = zeros(numel(funcList), 3);
curves = cell(numel(funcList), 1);

rng(2026, 'twister');

for i = 1:numel(funcList)
    funcNum = funcList(i);
    fobj = @(x) cec2017_eval(x, funcNum);
    [bestScore, bestPosition, curve] = MSDBO(popSize, maxIter, lb, ub, dim, fobj);

    results(i, :) = [funcNum, bestScore, norm(bestPosition)];
    curves{i} = curve;

    fprintf('F%02d  best = %.6e\n', funcNum, bestScore);
end

resultTable = array2table(results, 'VariableNames', {'Function', 'BestScore', 'BestPositionNorm'});
disp(resultTable);
writetable(resultTable, fullfile(rootDir, 'msdbo_cec2017_results.csv'));
save(fullfile(rootDir, 'msdbo_cec2017_curves.mat'), 'curves');
