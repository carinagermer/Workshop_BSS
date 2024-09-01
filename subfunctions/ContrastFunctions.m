function [G,DG] = ContrastFunctions(CF)
switch CF
    case "G0"
        % G = @(k) k.^4 / 4; % Kurtosis function
        G = @(x) 1/3*x.^3;
        DG = @(x) x.^2;
    case "G1"
        % G = @(k) k.^3 / 3; % Skewness function
        G = @(x) 1/2*x.^2;
        DG = @(x) x;
    case "G2"
        % G = @(k) log(cosh(k)); % Log-cosh function
        G = @(x) tanh(x);
        DG = @(x) 1-(tanh(x)).^2;
    case "G3"
        % G = @(k) exp(-k.^2 / 2); % Gaussian function
        G = @(x) -x.*exp(-x.^2/2);
        DG = @(x) exp(-x.^2/2).*(x.^2-1);
    case "G4"
        % G = @(k) log(k.^2 + 1); % Logarithm function
        G = @(x) 2.*x./(x.^2+1);
        DG = @(x) 2.*(1-x.^2)./(x.^2+1);
end
end

