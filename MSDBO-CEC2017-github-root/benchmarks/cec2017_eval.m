function f = cec2017_eval(x, funcNum)
%CEC2017_EVAL  Evaluate a CEC2017 function for one candidate solution.
%
%   f = cec2017_eval(x, funcNum)
%
%   If an official CEC2017 MATLAB/MEX evaluator named cec17_func is on the
%   MATLAB path, this wrapper calls it. Otherwise, it uses the lightweight
%   built-in benchmark functions in cec2017_basic.m.
%
%   For formal CEC2017 reporting, use the official CEC2017 benchmark files.

    if nargin ~= 2
        error('cec2017_eval requires x and funcNum.');
    end

    x = x(:)';

    if exist('cec17_func', 'file') == 2 || exist('cec17_func', 'file') == 3
        try
            f = cec17_func(x, funcNum);
        catch
            f = cec17_func(x(:), funcNum);
        end
        f = double(f(1));
        return;
    end

    f = cec2017_basic(x, funcNum);
end
