function [centre, radius] = stability_circles(s, inout,ax)
%STABILITY_CIRCLES Draw stability circles on a Smith Chart.
%   STABILITY_CIRCLES(S, INOUT) has the following parameters:
%   - S: a 2x2 matrix containing the s-parameters of the device, and
%   - INOUT: determines whether the input or output stability circles are
%     plotted.  The default is to plot the input stability circle.  The
%     options are:
%     -'i'   Plot the input stability cirle.
%     -'o'   Plot the output stability cirle.
%   If INOUT is not provided, the input stability circle will be plotted.
%   - ax: Axis to plot the circle(s) on, default: gca
%
%   Examples:
%   % The stability circles in Figure 3.3.6 in Gonzalez.
%   s500 =  [ (0.761*exp(-j*151/180*pi)) (0.025*exp( j* 31/180*pi)) ; ...
%             (11.84*exp( j*102/180*pi)) (0.429*exp(-j* 35/180*pi)) ]
%   s1000 = [ (0.770*exp(-j*166/180*pi)) (0.029*exp( j* 35/180*pi)) ; ...
%             ( 6.11*exp( j* 89/180*pi)) (0.365*exp(-j* 34/180*pi)) ]
%   smith
%   stability_circles(s500)
%   stability_circles(s1000)
%   figure
%   smith
%   stability_circles(s500, 'o')
%   stability_circles(s1000, 'o')
%   % The output stability circle in Figures 3.7.2 and 3.7.3 in Gonzalez.
%   s = [ 0.5*exp(-j*180/180*pi) 0.08*exp( j* 30/180*pi) ; ...
%         2.5*exp( j* 70/180*pi)  0.8*exp(-j*100/180*pi) ]
%   smith
%   stability_circles(s, 'o')
%
%   See also: SMITH, SMICIRC 
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.

% Check the number of inputs.
error(nargchk(1, 3, nargin));
% Check the s-parameter matrix has the correct size.
if ~isequal(size(s), [ 2 2 ])
  error('The s-parameter matrix must be a 2x2 matrix.');
end
% Check that the user has given a valid input.
if nargin > 1
  if (inout ~= 'i') && (inout ~= 'o')
    error('The second argument must be ''i'' or ''o''.');
  end
end

% Plot the input stability circle by default.
if nargin < 2
  inout = 'i';
end
if nargin < 3
  ax = gca;
end
% Add these circles to the current plot.
hold on

% Calculate the Rollett stability factor.
[ stable K delta ] = rollett_stability(s);

% Determine the centres and radii of the circles.
if ~stable
  if inout == 'i'
    centre = conj(s(1,1) - delta* ...
                  conj(s(2,2)))/(abs(s(1,1))^2 - abs(delta)^2);
    radius = abs(s(1,2)*s(2,1)/(abs(s(1,1))^2 - abs(delta)^2));
  else
    centre = conj(s(2,2) - delta* ...
                  conj(s(1,1)))/(abs(s(2,2))^2 - abs(delta)^2);
    radius = abs(s(1,2)*s(2,1)/(abs(s(2,2))^2 - abs(delta)^2));
  end
  smith_circles(centre, radius, 'r', ax)
end
if stable
    warning('MATLAB:stability_circles:Unconditionally_stable','Unconditionally Stable - No need to plot Stability Circle');
end

