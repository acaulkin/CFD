function PlotCp(Mesh, cp);
% function PlotCp(Mesh, cp);
% 
% This function plots the surface pressure contours on the
% body, given a mesh and cp at all elements
%
% INPUTS:
%   Mesh : mesh structure
%   cp   : nelem x 1 vector of pressure coefficients on all
%          elements
% OUTPUTS: a figure containing the cp distribution
%

hold on;

% we will plot on the body surface
BODY = 1;
B = Mesh.BC{BODY};
nedge = size(B,1);

for e=1:nedge, % loop over edges on body
  
  % info about this edge
  n1   = B(e,1); % first node
  n2   = B(e,2); % second node
  elem = B(e,3); % element
  
  % coordinates of nodes defining this edge
  x1 = Mesh.V(n1,:);
  x2 = Mesh.V(n2,:);
  
  % pressure coefficient on this edge
  cpe = cp(elem);
  
  % plot -cp vs. x (constant on edge)
  plot([x1(1), x2(1)], [-cpe, -cpe], 'b-', 'linewidth', 3);
  
end

set(gca, 'fontsize', 20);
grid on;
xlabel('x', 'fontsize', 24);
ylabel('-c_p', 'fontsize', 24);
