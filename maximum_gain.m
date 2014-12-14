function max_gain = maximum_gain(s)
%MAXIMUM_GAIN Calculate the maximum gain for a specified device.
%   MAX_GAIN = MAXIMUM_GAIN(S) has the following parameter:
%   - S: a 2x2 matrix containing the s-parameters of the device.
%   The output is:
%   - MAX_GAIN: The Maximum Available Gain (MAG) is computed if the device
%     is unconditionally stable, or the Maximum Stable Gain (MSG) is
%     computed if the device is not unconditionally stable.
%
%   Examples:
%   % The maximum gain in Example 3.6.1 in Gonzalez.
%   s = [ (0.641*exp(-j*171.3/180*pi)) (0.057*exp( j* 16.3/180*pi)) ; ...
%         (2.058*exp( j* 28.5/180*pi)) (0.572*exp(-j* 95.7/180*pi)) ]
%   10*log10(maximum_gain(s))  % Get the gain in decibels.
%   % The maximum gain in Example 3.7.3 in Gonzalez (plots in Figure
%   %  3.7.3).
%   s = [ (0.5*exp(-j*180/180*pi)) (0.08*exp( j* 30/180*pi)) ; ...
%         (2.5*exp( j* 70/180*pi)) ( 0.8*exp(-j*100/180*pi)) ]
%   10*log10(maximum_gain(s))  % Get the gain in decibels.
%
%   See also: ROLLETT_STABILITY, GA_CIRCLES, GP_CIRCLES
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.

% Check the number of inputs.
error(nargchk(1, 1, nargin));
% Check the s-parameter matrix has the correct size.
if ~isequal(size(s), [ 2 2 ])
  error('The s-parameter matrix must be a 2x2 matrix.');
end

% Calculate the Rollett stability factor.
[ stable K ] = rollett_stability(s);

% Calculate the maximum gain for each possible case.
if stable
    % The MAG because the device is unconditionally stable.
    max_gain = abs(s(2,1)/s(1,2))*(K - sqrt(K^2 - 1));
else
    % The MSG because the device is not unconditionally stable.
    max_gain = abs(s(2,1)/s(1,2));
end
