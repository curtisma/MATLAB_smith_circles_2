function values = gp_circles(s, values,ax)
%GP_CIRCLES Draw constant operating power gain circles on a Smith Chart.
%   GP_CIRCLES(S, VALUES) has the following parameters: 
%   - S: a 2x2 matrix containing the s-parameters of the device, and
%   - VALUES: the power values to be plotted in dB.
%   If VALUES is not provided, four circles are drawn from MAG or MSG in
%   steps of 1 dB.
%   - ax: Axis to plot the circle(s) on, default: gca
%
%   Example:
%   % The operating power gain circles in Figure 3.7.3 in Gonzalez.
%   s = [ (0.5*exp(-j*180/180*pi)) (0.08*exp( j* 30/180*pi)) ; ...
%         (2.5*exp( j* 70/180*pi)) ( 0.8*exp(-j*100/180*pi)) ]
%   smith
%   gp_circles(s, [ 30 20 15 10 0 ])
%
%   See also: SMITH, SMICIRC 
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.

% Check the number of inputs.
error(nargchk(1, 3, nargin));

% Calculate the values of the 4 default circles to be plotted if no gains
% are specified.
if nargin < 2
  values = 10.*log10(maximum_gain(s)) - [ 1e-4 1 2 3 ];
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
delta = s(1,1)*s(2,2) - s(1,2)*s(2,1);
K = (1 - abs(s(1,1))^2 - abs(s(2,2))^2 + abs(delta)^2)/ ...
    abs(2*s(1,2)*s(2,1));

% Determine the centres and radii of the circles.
gp = 10.^(values/10)/abs(s(2,1))^2;  % Divide by abs(s(2,1))^2 here to save
                                     % computations below.
centre = gp*conj(s(2,2) - delta*conj(s(1,1)))./ ...
         (1 + gp*(abs(s(2,2))^2 - abs(delta)^2));
radius = sqrt(1 - 2*K*abs(s(1,2)*s(2,1))*gp + (abs(s(1,2)*s(2,1))).^2*(gp).^2)./ ...
         abs((1 + gp*(abs(s(2,2))^2 - abs(delta)^2)));

% Plot the circles.
smith_circles(centre, abs(radius), 'm',256,ax)
