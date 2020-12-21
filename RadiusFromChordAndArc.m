function [R, delta] = RadiusFromChordAndArc(c, a)
% RadiusFromChordAndArc calculates the radius of a circle defined by a chord of
% length c and arc of length a. 
% 
% Equations to solve:
%   sin(delta) = c/a*delta
%   R = a/2/delta
% 
% [R, delta] = RadiusFromChordAndArc(c, a)
% 
% Inputs 
%   c           (numeric) [m]   Chord length
%   a           (numeric) [m]   Arc length
% Outputs
%   R           (numeric) [m]   Radius of the circle
%   delta       (numeric) [rad] Half-open angle of the circular arc
% 
% Shaun Fitzgerald

arguments
    c (:,:) double 
    a (:,:) double
end
assert(all(a > c), 'All input values of ''a'' must be greater than ''c''')

% The starting value is found by solving the second-degree Maclaurin polynomial
delta1 = sqrt(6*(1 - c./a));

cc = 0;
Err = 1;
while abs(Err) > 1e-14
    delta2 = delta1 - (sin(delta1) - c./a.*delta1)./(cos(delta1) - c./a);
    Err = sum(delta2 - delta1);
    delta1 = delta2;
    cc = cc + 1;
    if cc > 1000, error('Radius iteration not converged.'); end
end
delta = delta2;
R = a/2./delta; 