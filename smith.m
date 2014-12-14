function smith(z0, rcircles, xcircles, gcircles, bcircles)
%SMITH Draw a Smith chart.
%   SMITH(Z0, RCIRCLES, XCIRCLES, GCIRCLES, BCIRCLES) has the following
%   parameters:
%   - Z0: the normalising impedance (default 1 Ohm),
%   - RCIRCLES: the constant resistance circles to be drawn,
%   - XCIRCLES: the constant reactance circles to be drawn,
%   - GCIRCLES: the constant conductance circles to be drawn, and
%   - BCIRCLES: the constant susceptance circles to be drawn.
%   The circles will be normalised to Z0.  If no arguments are given, a
%   default Smith Chart will be plotted.
%
%   Examples:
%     smith
%     smith(50)
%     smith(1, [ 1/3 1 3 ], [ -3 -1 -1/3 1/3 1 3 ])
%     smith([], [], [], [ 1/3 1 3 ], [ -3 -1 -1/3 1/3 1 3 ])
%     smith(1, [ 1/3 1 3 ], [ -3 -1 -1/3 1/3 1 3 ], ...
%              [ 1/3 1 3 ], [ -3 -1 -1/3 1/3 1 3 ])
%     smith(1, [], [ 1/3 1 3 ], [], [ -3 -1 -1/3 ])
%
%   See also: smith_circles.m

% Check the number of inputs.
error(nargchk(0, 5, nargin));

% Set the default circles.
if nargin < 3
  xcircles = [];
end
if nargin < 4
  gcircles = [];
end
if nargin < 5
  bcircles = [];
end
if nargin == 0
  % The default Smith Chart.
  z0 = 1;
  rcircles = [ 1/3 1 3 ];
  xcircles = [ -3 -1 -1/3 1/3 1 3 ];
elseif nargin == 1
  % The default Smith Chart normalised to the specified z0.
  rcircles = z0*[ 1/3 1 3 ];
  xcircles = z0*[ -3 -1 -1/3 1/3 1 3 ];
end
if isempty(z0)
  z0 = 1;
end

% Normalise the circles for plotting
rvalues = rcircles/z0;
xvalues = xcircles/z0;
gvalues = gcircles*z0;
bvalues = bcircles*z0;

% A Smith Chart does not have the standard x and y axes that are used in
% a Matlab plot.
set(gca, 'XTick', [])        % Remove the x-axis labels.
set(gca, 'YTick', [])        % Remove the y-axis labels.
set(gca, 'XColor', 'w')      % Remove the x-axis.
set(gca, 'YColor', 'w')      % Remove the y-axis.
axis([ -1.1 1.1 -1.1 1.1 ])  % Make sure the labels are inside the axes.
% This will make the figure boundary square so that the Smith Chart is a
% circle rather than an elipse.
axis('square')
hold on

% Draw the constant conductance and susceptance circles.  Draw them first
% so that the constant resistance and reactance circles are plotted over
% them.
smith_circles(-[ gvalues./(gvalues + 1) (1 + j./bvalues) ], ...
              abs([ 1./(gvalues + 1) 1./bvalues ]), 'y--');

% Draw the horizontal line.
plot([ -1 1 ], [ 0 0 ], 'y');

% Add the r = 0 circle if not already present to draw the boundary of the
% Smith Chart.
if (size(rvalues, 2) == 0) || (min(abs(rvalues) ~= 0))
  rvalues = [ 0 rvalues ];
end

% Draw the constant resistance and reactance circles.
smith_circles([ rvalues./(rvalues + 1) (1 + j./xvalues) ], ...
              abs([ 1./(rvalues + 1) 1./xvalues ]), 'y');

% Label the circles.  Use the constant resistance and reactance circles
% unless they are not plotted.
if size(rcircles, 2) > 0
    
    % Label the constant resistance circles.
    xpos = (rcircles - z0)./(rcircles + z0);
    h = text(xpos', zeros(size(xpos')), num2str(rcircles', 2));
    set(h, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

    % Label zero and infinity.
    h = text(-1, 0, '0');
    set(h, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    h = text(1, 0, '\infty');
    set(h, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

else

    % Label the constant conductance circles.
    xpos = (1 - gcircles*z0)./(1 + gcircles*z0);
    h = text(xpos', zeros(size(xpos')), num2str(gcircles', 2));
    set(h, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

    % Label zero and infinity.
    h = text(-1, 0, '\infty');
    set(h, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    h = text(1, 0, '0');
    set(h, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

end

if size(xcircles, 2) > 0

  % Label the constant reactance circles.
  pos = (j*xcircles + z0)./(j*xcircles - z0);
  xpos = real(pos');
  ypos = imag(pos');
  h = text(xpos, ypos, [ num2str(xcircles', 2) 'j'*ones(size(xpos)) ]);
  set(h, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
  set(h(xpos > 0), 'HorizontalAlignment', 'left');
  set(h(xpos < 0), 'HorizontalAlignment', 'right');
  set(h(ypos > 0), 'VerticalAlignment', 'bottom');
  set(h(ypos < 0), 'VerticalAlignment', 'top');

else

  % Label the constant reactance circles.
  pos = (1 + j*z0*bcircles)./(1 - j*z0*bcircles);
  xpos = real(pos');
  ypos = imag(pos');
  h = text(xpos, ypos, [ num2str(bcircles', 2) 'j'*ones(size(xpos)) ]);
  set(h, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
  set(h(xpos > 0), 'HorizontalAlignment', 'left');
  set(h(xpos < 0), 'HorizontalAlignment', 'right');
  set(h(ypos > 0), 'VerticalAlignment', 'bottom');
  set(h(ypos < 0), 'VerticalAlignment', 'top');

end
