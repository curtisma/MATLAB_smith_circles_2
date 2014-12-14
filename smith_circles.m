function smith_circles(centre, radius, colour_linestyle, segments, ax)
%SMITH_CIRCLES Draw circles on a Smith Chart.
%   SMITH_CIRCLES(CENTRE, RADIUS, COLOUR_LINESTYLE, SEGMENTS) has the
%   following parameters:
%   - CENTRE: the centres of the circles (a complex number),
%   - RADIUS: the radii of the circles,
%   - COLOUR_LINESTYLE: one of the standard Matlab colour/linestyle
%     combinations that will be used to draw the circles, and
%   - SEGMENTS: the number of segments used to draw the circles.
%   - ax: Axis to plot the circle(s) on, default: gca
%   CENTRE and RADIUS must be row vectors with the same number of elements.
%   COLOUR_LINESTYLE and SEGMENTS are optional.  RADIUS must be positive and
%   purely real, and SEGMENTS must be an integer.
%
%   Examples:
%     smith_circle(0.5 + j*0.7, 1.5)
%     smith_circles([ 0 -0.5 0.5 -j*0.5 j*0.5 ], [ 1 0.5 0.5 0.5 0.5 ])
%     smith_circles(randn(1, 100)+j*randn(1, 100), randn(1, 100))
%     smith_circle(0.1 - j*0.3, 0.2, 'r--')
%     smith_circle(1.5 + j*1.1, 2.1, 'k:', 64)
%     smith_circles([ 0 -1 1 -j*1 j*1 ], [ 1 1 1 1 1 ], 'r-.', 2048)
%
%   See also: SMITH.M

% Check the number of inputs.
error(nargchk(2, 5, nargin));

% Set default values.
if nargin < 3
  % This is Matlab's default.
  colour_linestyle = 'b-';
end
if nargin < 4
  segments = 256;
end
if nargin < 5
  ax = gca;
end
% This values is used a number of times below, so let us calculate it
% once here rather than many times below.
ac = abs(centre);

% Find the circles that are entirely inside the Smith Chart.
inside = find(ac + radius <= 1);
% Find circles that intersect the Smith Chart boundary.
intersect = find((ac + radius > 1) & (abs(ac - radius) < 1));
% The other circles are all completely outside the Smith Chart or the
% Smith Chart is completely inside the circles - either way, they do not
% have to be plotted.

% The sizes of inside and intersect are used rather a lot, so they are
% stored in variables.
size_inside = numel(inside);
size_intersect = numel(intersect);

% Allocate memory for the points on the circles
z = zeros(segments + 1, size_inside + size_intersect);

% Circles that are entirely inside the Smith Chart
if size_inside > 0

  % The angles from the centres of the circles to the points to be plotted.  I
  % have to add 1 to the number of segments because there is one more point
  % than there are segments.
  theta = linspace(-pi, pi, segments + 1)';

  % Find the points on the circles.
  z(:, (1:size_inside)) = ones(size(theta))*centre(inside) + ...
      exp(j*theta)*radius(inside);

end

% Circles that intersect the Smith Chart boundary.
if size_intersect > 0

  % Find the angles to the intercept points.
  % The angle from the centre to the intercept points relative to the angle to
  % the centre is calculated.  This calculation is based on the well-known
  % relationship:
  % c^2 = a^2 + b^2 - 2*a*b*cos(theta)
  % Note that there is a -1 factor implicit in the next line which
  % ensures that the correct angle is used.
  theta_p = acos((1 - radius(intersect).^2 - ac(intersect).^2)./ ...
                 (2.*radius(intersect).*ac(intersect)));
  theta_1 = angle(centre(intersect)) - theta_p;
  theta_2 = angle(centre(intersect)) + theta_p;
  % I need to change the order of these angles to get the part of the circle
  % that is inside the Smith Chart boundary.
  theta_1 = theta_1 + 2*pi;

  % The angles from the centres of the circles to the points to be plotted.  I
  % have to add 1 to the number of segments because there is one more point
  % than there are segments.
  theta = ones(segments + 1, 1)*theta_1 + ...
          linspace(0, 1, segments + 1)'*(theta_2 - theta_1);

  % Find the points on the circles.
  z(:, ((size_inside + 1):(size_inside + size_intersect))) = ...
      ones(segments + 1, 1)*centre(intersect) + ...
      exp(j*theta).*(ones(segments + 1, 1)*radius(intersect));

end

% Finally, plot all the circles.
plot(ax, z, colour_linestyle);
