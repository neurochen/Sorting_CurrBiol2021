function y = betad(x,a,b);
%BETAd	Beta distribution function.
%	BETAd(x,a,b) 

%	Gordon Smyth, gks@maths.uq.edu.au, University of Queensland
%	1 August 1998

% Initialize Y to zero.
y = zeros(size(x));

% Return NaN for parameter values outside their respective limits.
k1 = find(a <= 0 | b <= 0 | x < 0 | x > 1);
if any(k1)
    tmp = NaN;
    y(k1) = tmp(ones(size(k1))); 
end

% Return Inf for x = 0 and a < 1 or x = 1 and b < 1.
% Required for non-IEEE machines.
k2 = find((x == 0 & a < 1) | (x == 1 & b < 1));
if any(k2)
    tmp = Inf;
    y(k2) = tmp(ones(size(k2))); 
end

% Return the beta density function for valid parameters.
k = find(~(a <= 0 | b <= 0 | x <= 0 | x >= 1));
if any(k)
%    y(k) = x(k) .^ (a(k) - 1) .* (1 - x(k)) .^ (b(k) - 1) ./ beta(a(k),b(k));
     tmp(k) = (a - 1).*log(x(k)) + (b - 1).*log((1 - x(k))) - betaln(a,b);
     y(k) = exp(tmp(k));
end

