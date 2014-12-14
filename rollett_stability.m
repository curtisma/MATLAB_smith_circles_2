function [ stable, K, delta ] = rollett_stability(s)
%ROLLETT_STABILITY Evaluate the Rollett stability criteria.
%   [ STABLE, K, DELTA ] = ROLLETT_STABILITY(S) has the following
%   parameter:
%   - S: a 2x2 matrix containing the s-parameters of the device.
%   The outputs are:
%   - STABLE: 1 if the device is unconditionally stable and 0 if not (see
%     the note below),
%   - K: the value of the Rollett Stability Factor, and
%   - DELTA: the value of delta used to caculate K and test for
%     unconditional stability.
%
%   NOTE: Only considering K and DELTA (as is almost universally done -
%   e.g. Gonzalez) can lead to the incorrect conclusion being reached about
%   whether a device is unconditionally stable or not.  Such problems are
%   rare, but can occur in practice.
%
%   Examples:
%   % Check the stability of the device in Example 3.3.1 in Gonzalez
%   % (plotted in  Figure 3.3.6).
%   s500 =  [ (0.761*exp(-j*151/180*pi)) (0.025*exp( j* 31/180*pi)) ; ...
%             (11.84*exp( j*102/180*pi)) (0.429*exp(-j* 35/180*pi)) ]
%   s1000 = [ (0.770*exp(-j*166/180*pi)) (0.029*exp( j* 35/180*pi)) ; ...
%             ( 6.11*exp( j* 89/180*pi)) (0.365*exp(-j* 34/180*pi)) ]
%   s2000 = [ (0.760*exp(-j*174/180*pi)) (0.040*exp( j* 44/180*pi)) ; ...
%             ( 3.06*exp( j* 74/180*pi)) (0.364*exp(-j* 43/180*pi)) ]
%   s4000 = [ (0.756*exp(-j*179/180*pi)) (0.064*exp( j* 48/180*pi)) ; ...
%             ( 1.53*exp( j* 53/180*pi)) (0.423*exp(-j* 66/180*pi)) ]
%   [ stable500  K500  delta500 ]  = rollett_stability(s500)
%   [ stable1000 K1000 delta1000 ] = rollett_stability(s1000)
%   [ stable2000 K2000 delta2000 ] = rollett_stability(s2000)
%   [ stable4000 K4000 delta4000 ] = rollett_stability(s4000)
%
%   See also: STABILITY_CIRCLES, GA_CIRCLES, GP_CIRCLES
%
%   Reference:
%   G. Gonzalez, "Microwave Transistor Amplifiers - Analysis and Design,"
%   2nd ed., Prentice Hall, 1997.
%   R. W. Jackson, "Rollett Proviso in the Stability of Linear Microwave
%   Circuits - A Tutorial," IEEE Transactions on Microwave Theory and
%   Techniques, Vol. 54, No. 3, March 2006, pp. 993-1000.

% Check the number of inputs.
error(nargchk(1, 1, nargin));
% Check the s-parameter matrix has the correct size.
if ~isequal(size(s), [ 2 2 ])
  error('The s-parameter matrix must be a 2x2 matrix.');
end

% Calculate the Rollett stability factor.
delta = s(1,1)*s(2,2) - s(1,2)*s(2,1);
K = (1 - abs(s(1,1))^2 - abs(s(2,2))^2 + abs(delta)^2)/ ...
    abs(2*s(1,2)*s(2,1));

% Check for unconditional stability.
stable = 0;
if (K > 1) && (abs(delta) < 1)
    stable = 1;
end
