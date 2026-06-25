function f = cec2017_basic(x, funcNum)
%CEC2017_BASIC  Lightweight CEC2017-style benchmark functions.
%
%   This file is included so that the repository can run without external
%   benchmark files. It is not a replacement for the official shifted and
%   rotated CEC2017 implementation.

    x = x(:)';
    D = numel(x);

    switch funcNum
        case 1
            f = elliptic(x);
        case 2
            f = bentCigar(x);
        case 3
            f = discus(x);
        case 4
            f = rosenbrock(x);
        case 5
            f = ackley(x);
        case 6
            f = weierstrass(x);
        case 7
            f = griewank(x);
        case 8
            f = rastrigin(x);
        case 9
            f = schwefel(x);
        case 10
            f = katsuura(x);
        case 11
            f = happyCat(x);
        case 12
            f = hgbat(x);
        case 13
            f = expandedGriewankRosenbrock(x);
        case 14
            f = expandedScafferF6(x);
        case 15
            f = levy(x);
        case 16
            f = zakharov(x);
        case 17
            f = sumDifferentPowers(x);
        case 18
            f = dixonPrice(x);
        case 19
            f = griewank(x) + 0.1 * rastrigin(x);
        case 20
            f = ackley(x) + 0.1 * elliptic(x);
        case 21
            f = hybridFunction(x, {@rastrigin, @ackley, @schwefel});
        case 22
            f = hybridFunction(x, {@griewank, @rosenbrock, @weierstrass});
        case 23
            f = hybridFunction(x, {@happyCat, @hgbat, @rastrigin});
        case 24
            f = hybridFunction(x, {@levy, @ackley, @griewank});
        case 25
            f = hybridFunction(x, {@katsuura, @rastrigin, @schwefel});
        case 26
            f = compositionFunction(x, {@rastrigin, @ackley, @griewank}, [0, 100, 200]);
        case 27
            f = compositionFunction(x, {@rosenbrock, @schwefel, @elliptic}, [0, 100, 200]);
        case 28
            f = compositionFunction(x, {@weierstrass, @griewank, @ackley}, [0, 100, 200]);
        case 29
            f = compositionFunction(x, {@happyCat, @hgbat, @rastrigin}, [0, 100, 200]);
        case 30
            f = compositionFunction(x, {@katsuura, @levy, @schwefel}, [0, 100, 200]);
        otherwise
            error('funcNum must be an integer from 1 to 30.');
    end

    if ~isfinite(f)
        f = realmax;
    end
end

function f = elliptic(x)
    D = numel(x);
    if D == 1
        f = x.^2;
    else
        weights = (1e6) .^ ((0:D-1) / (D-1));
        f = sum(weights .* x.^2);
    end
end

function f = bentCigar(x)
    f = x(1)^2 + 1e6 * sum(x(2:end).^2);
end

function f = discus(x)
    f = 1e6 * x(1)^2 + sum(x(2:end).^2);
end

function f = rosenbrock(x)
    z = x + 1;
    f = sum(100 * (z(2:end) - z(1:end-1).^2).^2 + (z(1:end-1) - 1).^2);
end

function f = ackley(x)
    D = numel(x);
    f = -20 * exp(-0.2 * sqrt(sum(x.^2) / D)) - exp(sum(cos(2 * pi * x)) / D) + 20 + exp(1);
end

function f = weierstrass(x)
    a = 0.5;
    b = 3;
    kMax = 20;
    D = numel(x);
    f = 0;
    for i = 1:D
        for k = 0:kMax
            f = f + a^k * cos(2 * pi * b^k * (x(i) + 0.5));
        end
    end
    c = 0;
    for k = 0:kMax
        c = c + a^k * cos(2 * pi * b^k * 0.5);
    end
    f = f - D * c;
end

function f = griewank(x)
    D = numel(x);
    f = 1 + sum(x.^2) / 4000 - prod(cos(x ./ sqrt(1:D)));
end

function f = rastrigin(x)
    D = numel(x);
    f = sum(x.^2 - 10 * cos(2 * pi * x) + 10) + 10 * D;
end

function f = schwefel(x)
    D = numel(x);
    z = x + 420.9687462275036;
    f = 418.9829 * D - sum(z .* sin(sqrt(abs(z))));
end

function f = katsuura(x)
    D = numel(x);
    prodTerm = 1;
    for i = 1:D
        sumTerm = 0;
        for j = 1:32
            sumTerm = sumTerm + abs(2^j * x(i) - round(2^j * x(i))) / 2^j;
        end
        prodTerm = prodTerm * (1 + (i * sumTerm))^(10 / D^1.2);
    end
    f = (10 / D^2) * prodTerm - 10 / D^2;
end

function f = happyCat(x)
    D = numel(x);
    alpha = 1 / 8;
    r2 = sum(x.^2);
    f = (abs(r2 - D))^(2 * alpha) + (0.5 * r2 + sum(x)) / D + 0.5;
end

function f = hgbat(x)
    D = numel(x);
    r2 = sum(x.^2);
    sx = sum(x);
    f = sqrt(abs(r2^2 - sx^2)) + (0.5 * r2 + sx) / D + 0.5;
end

function f = expandedGriewankRosenbrock(x)
    D = numel(x);
    z = x + 1;
    f = 0;
    for i = 1:D
        j = i + 1;
        if j > D
            j = 1;
        end
        temp = 100 * (z(i)^2 - z(j))^2 + (z(i) - 1)^2;
        f = f + temp^2 / 4000 - cos(temp) + 1;
    end
end

function f = expandedScafferF6(x)
    D = numel(x);
    f = 0;
    for i = 1:D
        j = i + 1;
        if j > D
            j = 1;
        end
        s = x(i)^2 + x(j)^2;
        f = f + 0.5 + (sin(sqrt(s))^2 - 0.5) / (1 + 0.001 * s)^2;
    end
end

function f = levy(x)
    w = 1 + (x - 1) / 4;
    f = sin(pi * w(1))^2 + sum((w(1:end-1) - 1).^2 .* (1 + 10 * sin(pi * w(1:end-1) + 1).^2)) + ...
        (w(end) - 1)^2 * (1 + sin(2 * pi * w(end))^2);
end

function f = zakharov(x)
    D = numel(x);
    i = 1:D;
    s = sum(0.5 * i .* x);
    f = sum(x.^2) + s^2 + s^4;
end

function f = sumDifferentPowers(x)
    D = numel(x);
    f = sum(abs(x) .^ (2 + 4 * ((1:D) - 1) / max(D - 1, 1)));
end

function f = dixonPrice(x)
    D = numel(x);
    if D == 1
        f = (x(1) - 1)^2;
    else
        i = 2:D;
        f = (x(1) - 1)^2 + sum(i .* (2 * x(2:end).^2 - x(1:end-1)).^2);
    end
end

function f = hybridFunction(x, funcs)
    D = numel(x);
    parts = splitIndex(D, numel(funcs));
    f = 0;
    for k = 1:numel(funcs)
        idx = parts{k};
        if ~isempty(idx)
            f = f + funcs{k}(x(idx));
        end
    end
end

function f = compositionFunction(x, funcs, bias)
    m = numel(funcs);
    values = zeros(1, m);
    weights = zeros(1, m);
    for k = 1:m
        shift = (k - 1) * 5;
        z = x - shift;
        values(k) = funcs{k}(z) + bias(k);
        weights(k) = 1 / (sqrt(sum(z.^2)) + 1e-12);
    end
    weights = weights / sum(weights);
    f = sum(weights .* values);
end

function parts = splitIndex(D, nParts)
    edges = round(linspace(0, D, nParts + 1));
    parts = cell(1, nParts);
    for k = 1:nParts
        parts{k} = (edges(k) + 1):edges(k + 1);
    end
end
