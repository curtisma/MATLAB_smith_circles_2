function noise_circles(fmin, refo, rn, values, ax)
%NOISE_CIRCLES Draw constant noise figure circles on a Smith Chart.
%   NOISE_CIRCLES(FMIN, REFO, RN, VALUES) has the following parameters:
%   - FMIN: the minimum noise figure in decibels (dB),
%   - GAMMA_0: the optimum  input  reflection  coefficient,
%   - RN: the normalised  noise resistance, and 
%   - VALUES: the noise figures for which circles are to be plotted in
%   decibels (dB).
%   If VALUES is not provided, four circles are drawn from FMIN in steps of
%   0.5 dB.
%   - ax: Axis to plot the circle(s) on, default: gca
%
%   Example:
%   % The available power gain circles in Figure 4.3.4 in Gonzalez.
%   fmin = 2.5
%   gamma_opt = 0.475*exp(j*166/180*pi)
%   rn = 3.5/50
%   smith
%   noise_circles(fmin, gamma_opt, rn, (2.5:0.1:3) + 1e-5)
%   % The small value is added to make the smallest circle visible.
%
%   See also: SMITH, SMITH_CIRLES 
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.

% Check the number of inputs.
error(nargchk(3, 5, nargin));

% Check for the minimum number of inputs.
if nargin < 3
  error('There must be at least 3 inputs.');
end
% The values for the 4 default circles that will be plotted if no noise
% figures are specified.
if nargin < 4
  values = fmin + [1e-4 0.5 1 1.5];
end
if nargin < 5
  ax = gca;
end
% Add these circles to the current plot.
hold on

% Determine the centres and radii of the circles.
N = (10.^(values/10) - 10^(fmin/10))/(4*rn)*abs(1 + refo)^2;
centre = refo./(1 + N);
radius = abs(sqrt(N.^2 + N*(1 - abs(refo)^2))./(1 + N));

% Plot the circles.
smith_circles(centre, radius, 'b',256,ax)
