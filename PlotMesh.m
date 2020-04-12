function PlotMesh(Mesh);
% function PlotMesh(Mesh);
%
% PURPOSE: Plots the mesh
%
% INPUTS:
%   Mesh    : Mesh to plot
%
% The input Mesh is a structure with the following fields:
%
% Mesh.V   = nNode x dim array of node coordinates.  For each node i,
%            V(i,:) gives the x and y location of the node.
%
% Mesh.E2N = nElem x 3 array of nodes for each element.  For each
%            element k, E2N(k,:) is a list of the 3 node numbers that
%            form the element.
%
% Mesh.BC  = Cell array of 6 matrices corresponding to the 6
%            boundaries.  The names of these boundaries are given in
%            Mesh.BCTitle.  Each matrix contains nEdge rows and 3
%            columns.  nEdge is the number of edges on the boundary.
%            The 3 columns contain [n1, n2, elem], where n1 and n2 are
%            the node numbers that define the edge, and elem is the
%            index of the adjacent element.  Note, n1 and n2 are order
%            counter-clockwise from the point of view of the adjacent
%            element.
%
% Mesh.BCTitle = Cell array of strings identifying the boundary names.
%

%--------------------------------------
% plot using pdeplot
t = Mesh.E2N';
t = [t; ones(1,size(t,2))];
h=pdeplot(Mesh.V', [],t,'mesh','on');
set(h,'color', [.5, .5, .5]);

%-------------------------------------
% set zoom
axis equal;
axis off;
xyrange = [-2.2, 3.2, -0.05, 1.05];
axis(xyrange);
AR = (xyrange(2)-xyrange(1))/(xyrange(4)-xyrange(3));
DY = 400; % plot size in y (not relevant for width/height of export)
set(gcf, 'position', [100, 100, DY*AR, DY]);
set(gca, 'position', [0,0,1,1]);
set(gcf, 'paperposition', [0.25, 2.5, 8, 8/AR]);
