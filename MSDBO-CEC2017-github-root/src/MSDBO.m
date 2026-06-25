function [fMin, bestX, Convergence_curve] = MSDBO(N, iter, lb, ub, dim, fobj)
%MSDBO  Minimal MATLAB implementation of MSDBO.
%
%   [fMin, bestX, Convergence_curve] = MSDBO(N, iter, lb, ub, dim, fobj)
%
%   Inputs
%   -----
%   N      : population size. Use N >= 20.
%   iter   : maximum number of iterations.
%   lb     : scalar or 1-by-dim lower bound.
%   ub     : scalar or 1-by-dim upper bound.
%   dim    : problem dimension.
%   fobj   : objective function handle. It receives a row vector x.
%
%   Outputs
%   -------
%   fMin              : best objective value found.
%   bestX             : best solution found.
%   Convergence_curve : best objective value at initialization and each iteration.

    if nargin ~= 6
        error('MSDBO requires six inputs: N, iter, lb, ub, dim, fobj.');
    end
    if N < 20
        error('MSDBO uses fixed update groups and requires N >= 20.');
    end
    if iter < 1
        error('iter must be a positive integer.');
    end

    lb = expandBound(lb, dim);
    ub = expandBound(ub, dim);
    if any(lb >= ub)
        error('Each lower bound must be smaller than the corresponding upper bound.');
    end

    x = rand(N, dim) .* (ub - lb) + lb;
    fit = zeros(1, N);
    for i = 1:N
        fit(i) = fobj(x(i, :));
    end

    pX = x;
    pFit = fit;
    previousPX = pX;

    [fMin, bestIndex] = min(fit);
    bestX = x(bestIndex, :);
    initialBest = fMin;
    Convergence_curve = zeros(1, iter + 1);
    Convergence_curve(1) = initialBest;

    for t = 1:iter
        pPercentLow = 0.10;
        pPercentHigh = 0.30;
        pPercent = pPercentLow + (pPercentHigh - pPercentLow) * t / iter;
        pNum = round(N * pPercent);

        [~, sortIndex] = sort(fit);
        bestX2 = pX(sortIndex(2), :);

        [~, worstIndex] = max(fit);
        worstX = x(worstIndex, :);
        r2 = rand;

        for i = 1:pNum
            if r2 < 0.9
                direction = 1;
                if rand <= 0.1
                    direction = -1;
                end
                x(i, :) = pX(i, :) + 0.3 * abs(pX(i, :) - worstX) + direction * 0.1 * previousPX(i, :);
            else
                angle = randperm(180, 1);
                if angle == 90 || angle == 180
                    x(i, :) = pX(i, :);
                else
                    theta = angle * pi / 180;
                    x(i, :) = pX(i, :) + tan(theta) .* abs(pX(i, :) - previousPX(i, :));
                end
            end
            x(i, :) = newBound(x(i, :), lb, ub, bestX, bestX2);
            fit(i) = fobj(x(i, :));
        end

        [~, currentBestIndex] = min(fit);
        currentBestX = x(currentBestIndex, :);

        R = 1 - t / iter;
        Xnew1 = applyBounds(currentBestX .* (1 - R), lb, ub);
        Xnew2 = applyBounds(currentBestX .* (1 + R), lb, ub);
        Xnew11 = applyBounds(bestX .* (1 - R), lb, ub);
        Xnew22 = applyBounds(bestX .* (1 + R), lb, ub);

        for i = (pNum + 1):12
            x(i, :) = currentBestX + rand(1, dim) .* (pX(i, :) - Xnew1) + rand(1, dim) .* (pX(i, :) - Xnew2);
            x(i, :) = newBound(x(i, :), lb, ub, bestX, bestX2);
            fit(i) = fobj(x(i, :));
        end

        for i = 13:19
            if rand < 0.5
                x(i, :) = pX(i, :) + randn .* (pX(i, :) - Xnew11) + rand(1, dim) .* (pX(i, :) - Xnew22);
            else
                wMin = 1;
                wMax = 3.0;
                w = wMin + rand * (wMax - wMin);
                w = abs(((2 * rand) - (rand + rand)) / w^2);
                x(i, :) = rand * pX(i, :) + w * (bestX - pX(i, :));
            end
            x(i, :) = newBound(x(i, :), lb, ub, bestX, bestX2);
            fit(i) = fobj(x(i, :));
        end

        for j = 20:N
            x(j, :) = bestX + randn(1, dim) .* (abs(pX(j, :) - currentBestX) + abs(pX(j, :) - bestX)) ./ 2;
            x(j, :) = newBound(x(j, :), lb, ub, bestX, bestX2);
            fit(j) = fobj(x(j, :));
        end

        previousPX = pX;
        for i = 1:N
            if fit(i) < pFit(i)
                pFit(i) = fit(i);
                pX(i, :) = x(i, :);
            end
            if pFit(i) < fMin
                fMin = pFit(i);
                bestX = pX(i, :);
            end
        end

        Convergence_curve(t + 1) = fMin;
    end
end

function bound = expandBound(bound, dim)
    if isscalar(bound)
        bound = repmat(bound, 1, dim);
    else
        bound = bound(:)';
        if numel(bound) ~= dim
            error('Bound vector length must match dim.');
        end
    end
end

function x = applyBounds(x, lb, ub)
    x = min(max(x, lb), ub);
end

function x = newBound(x, lb, ub, best1, best2)
    lowMask = x < lb;
    x(lowMask) = rand * best1(lowMask) + rand * (best1(lowMask) - best2(lowMask));

    highMask = x > ub;
    x(highMask) = rand * best2(highMask) + rand * (best1(highMask) - best2(highMask));

    x = applyBounds(x, lb, ub);
end
