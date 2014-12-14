function values = ga_circles(s, values, ax)
%GA_CIRCLES Draw constant available power gain circles on a Smith Chart.
%   GA_CIRCLES(S, VALUES) has the following parameters:
%   - S: a 2x2 matrix containing the s-parameters of the device, and
%   - VALUES: the power values to  be plotted in dB.
%   If VALUES is not provided, four circles are drawn from MAG or MSG in
%   steps of 1 dB.
%   - ax: Axis to plot the circle(s) on, default: gca
%
%   Example:
%   % The available power gain circles in Figure 4.3.4 in Gonzalez.
%   s = [ (0.552*exp( j*169/180*pi)) (0.049*exp( j* 23/180*pi)) ; ...
%         (1.681*exp( j* 26/180*pi)) (0.839*exp(-j* 67/180*pi)) ]
%   smith
%   ga_circles(s, 13.7:-1:6.7)
%
%   See also: SMITH, SMITH_CIRLES 
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.

% Check the number of inputs.
error(nargchk(1, 3, nargin));

% Calculate the values of the 4 default circles to be plotted if no gains
% are specified.
if nargin < 2
  values = 10*log10(maximum_gain(s)) - [ 1e-4 1 2 3 ];
end
if nargin < 3
  ax = gca;
end

% Check the s-parameter matrix has the correct size.
if ~isequal(size(s), [ 2 2 ])
  error('The s-parameter matrix must be a 2x2 matrix.');
end

% Add these circles to the current plot.
hold on

% Calculate the Rollett stability factor.
[ stable K delta ] = rollett_stability(s);

% Determine the centres and radii of the circles.
ga = 10.^(values/10)/abs(s(2,1))^2;  % Divide by abs(s(2,1))^2 here to save
                                     % computations below.
centre = ga*conj(s(1,1) - delta*conj(s(2,2)))./ ...
         (1 + ga*(abs(s(1,1))^2 - abs(delta)^2));
radius = sqrt(1 - 2*K*abs(s(1,2)*s(2,1))*ga + (abs(s(1,2)*s(2,1))*ga).^2)./ ...
         (1 + ga*(abs(s(1,1))^2 - abs(delta)^2));

% Plot the Ga circles.
smith_circles(centre, abs(radius), 'm', 256, ax)
