%Andrew Caulkins
%Project 3 Script
%___________________________________________
clc
clear
close all

syms c eta

%Loading in Mesh
Mesh = load('Mesh2.mat');


%Making Mesh more "usable" by converting from structure to a matrix:
Mesh = cell2mat(struct2cell(Mesh));

% 
% %Plotting Mesh in Order to Visualize System:
figure(1)
PlotMesh(Mesh)
title('Finite-Element System')

%Length of System:
L = 5.0;


%% Determining w(x) at each node in a given mesh
%it is first important to know what the basis functions in the mesh are,
%So I will compute those now:


%Determine the x and y coordinates of each node in the mesh, and the
%element associated with each:

Elements = Mesh.E2N;
Node_Locs = Mesh.V;

%The idea is to set up a jacobian matrix for each element, and solve for
%the basis function for each element

%Setting up column vectors for each node inside element:
Element_1 = Elements(:,1);
Element_2 = Elements(:,2);
Element_3 = Elements(:,3);

%Determining Element Node Locations:
Element_Loc_1 = Node_Locs(Element_1,:);
Element_Loc_2 = Node_Locs(Element_2,:);
Element_Loc_3 = Node_Locs(Element_3,:);

Element_Loc_1X = Element_Loc_1(:,1);
Element_Loc_2X = Element_Loc_2(:,1);
Element_Loc_3X = Element_Loc_3(:,1);

Element_Loc_1Y = Element_Loc_1(:,2);
Element_Loc_2Y = Element_Loc_2(:,2);
Element_Loc_3Y = Element_Loc_3(:,2);



%Putting Into Single Array Containing the x and y positions of each node
%composing an element:
Element_Locations = [Element_Loc_1,Element_Loc_2,Element_Loc_3];

%Number of elements:
N = size(Element_Locations,1);

%Initializing Basis Functions:
phi_1 = 1 - c - eta;
phi_2 = c;
phi_3 = eta;
phi = [phi_1 phi_2 phi_3];


%Preallocating Sparse Matrix:
A_global = sparse(size(Node_Locs,1),size(Node_Locs,1));

%Constructing Local A Matrix:
for i=1:size(Element_Locations,1)
    I = Mesh.E2N(i,:);
    %First Calculating Local Jacobian:
    J = Element_Locations(i,:);
    J_bar = [J(3)-J(1) , J(5)-J(1) ;
        J(4)-J(2) , J(6)-J(2)];
    J_det = det(J_bar);
   
   
    Area = 0.5*J_det;
    derivs = [-1 -1 ; 1 0 ; 0 1];
    Gradient = derivs*inv(J_bar);
    
    
   %Calculating LOCAL A matrix:
    for O=1:3
        for M=1:3
           A_k(O,M) = -dot(Gradient(O,:),Gradient(M,:))*Area;
        end 
    end
    
   %Assembling GLOBAL A matrix:
        A_global(I,I) = A_global(I,I) + A_k;
        F(I) = 0;
end

    %Applying Inflow Boundary Conditions:
    for i=1:size(Mesh.BC{5},1)
        I = Mesh.BC{5}(i,1);
        A_global(I,:) = 0;
        A_global(I,I) = diag(1);
        F(I) = 0;
    end
    I = Mesh.BC{5}(end,2);
    A_global(I,:) = 0;
    A_global(I,I) = diag(1);
    F(I) = 0;
        
    
    %Applying Outflow Boundary Conditions:
    for i=1:size(Mesh.BC{3},1)
        I = Mesh.BC{3}(i,1);
        A_global(I,:) = 0;
        A_global(I,I) = diag(1);
        F(I) = L;
    end
    I = Mesh.BC{3}(end,2);
    A_global(I,:) = 0;
    A_global(I,I) = diag(1);
    F(I) = L;

            
P = 1;
w = -A_global\F';
%Calculating w(x) from known values:
for i=1:size(Element_Locations,1)
    I = Mesh.E2N(i,:);
    %First Calculating Local Jacobian:
    J = Element_Locations(i,:);
    J_bar = [J(3)-J(1) , J(5)-J(1) ;
        J(4)-J(2) , J(6)-J(2)];
    J_det = det(J_bar);
   
   
    Area = 0.5*J_det;
    derivs = [-1 -1 ; 1 0 ; 0 1];
    Gradient = derivs*inv(J_bar);
    
    
    
    %Looping to Determine gradient of w:
    for j=1:3
        delta_w(j,:) = w(I(j))*Gradient(j,:);
    end
    
    Gradient_w(i,:) = sum(delta_w);
    mag(i) = norm(Gradient_w(i,:));
end
    
    %Plotting
    figure(2)
    pdeplot(Mesh.V', [], [Mesh.E2N'; ones(1,size(Element_Locations,1))], 'XYData', mag)
    xlabel('x-Distance [m]')
    ylabel('y-Distance [m]')
    title('Mesh Velocity Contours')
  

%% Post Processing 

%Body:
Body = Mesh.BC(1);
Body = cell2mat(Body);
Body_3 = Body(:,3);
%Node Coordinates:
T = Mesh.E2N(Body_3,:);

%Plotting Body Boundary:
figure(3)
t = T';
t = [t; ones(1,size(t,2))];
h=pdeplot(Mesh.V', [],t,'mesh','on'); 
set(h,'color', [.5, .5, .5]);
title('Body Inside Flow')
xlabel('Length Along Body [m]')
ylabel('y- Distance [m]')
hold on
    
Edges = zeros(size(T,1)-1,1);
%Determining Edge Node Numbers:
for i=1:size(Edges,1)
    Edges(i) = min(intersect(T(i,:),T(i+1,:)));
end



%Determining X,Y Location of Each Found Node:
Edges_Loc = Mesh.V(Edges,:);
Edges_LocX = Edges_Loc(:,1);
Edges_LocY = Edges_Loc(:,2);

%Determining Normals and Tangents:
x = Edges_Loc(:,1); y = Edges_Loc(:,2); % Node x and y Locations
xm = 0.5*(x(1:end-1) + x(2:end)); % interval midpoints: x
ym = 0.5*(y(1:end-1) + y(2:end)); % interval midpoints: y
% xm = [1;xm;1]; ym = [0;ym;0]; % include trailing edge as endpoint
dx = xm(2:end) - xm(1:end-1); % interval delta x
dy = ym(2:end) - ym(1:end-1); % interval delta y
ds = sqrt(dx.^2 + dy.^2); % interval length
t = [dx./ds,dy./ds]; % tangent vector
n = [dy./ds,-dx./ds]; % normal vector

xm(end) = [];
ym(end) = []; 


% %Plotting:
scale = 0.25;
quiver(xm,ym,n(:,1),n(:,2),scale)
quiver(xm,ym,t(:,1),t(:,2),scale)
legend('Body Boundary','Normal Vectors','Tangent Vectors')



%Explain how I will determine the pressure Coefficient on each edge point

%Determining Velocity Values at each Body Node:
Element_Numbers = Body(:,3);
Body_mags = Gradient_w(Element_Numbers,:);

%Calculating Coefficient of Pressure:
c_pe = 1-(Body_mags(:,1).^2 + Body_mags(:,2).^2);

%Plotting cp Values:
figure(4)
X = linspace(0,1,size(c_pe,1));
plot(X,c_pe)
title('Coefficient of Pressure')
axis ij  %Reversing Axis
       

js = zeros(size(ds,1),2);
js(:,2) = 1;
Dists = ds .* js;


%Determining Sectional Lift:
delta_p = .5*1.225 *c_pe;
delta_p(end) = [];
delta_p(end) = [];
delta_p(end) = [];
InnerX = -delta_p.*n(:,1);
InnerY = -delta_p.*n(:,2);

Inner = [InnerX,InnerY];
for i=1:size(Inner,1)
L_e(i) = -dot(Inner(i,:),Dists(i,:));
end



figure(5)
[hAx,hLine1,hLine2] = plotyy(X,c_pe,x,y);
legend('Coefficient of Pressure','Body Boundary')
title('Cp Overlay')
xlabel('Length Along Body [m]')
ylabel(hAx(1),'Coefficient of Pressure [unitless]') % left y-axis 
ylabel(hAx(2),'Distance [m]') % right y-axis
axis ij  %Reversing Axis



%Computing c_l:
c_l = trapz(L_e)/((1/2)*1.225)



%% Determining Coefficient of Lift by Calculating Coefficient of Pressure along top wall and Symmetry Boundaries:
%Determining Velocity Values at each Body Node:

%Top Wall:
Top = Mesh.BC(4);
Top = cell2mat(Top);

%Computing Magnitude of Velocity along top wall:
Element_Numbers_Top = Top(:,3);
Top_mags = Gradient_w(Element_Numbers_Top,:);

%Calculating Coefficient of Pressure Along top wall:
c_p_Top = 1-(Top_mags(:,1).^2 + Top_mags(:,2).^2);


%Symmetry Right:
Right = Mesh.BC(2);
Right = cell2mat(Right);

%Computing Magnitude of Velocity along Right Symmetry:
Element_Numbers_Right = Right(:,3);
Right_mags = Gradient_w(Element_Numbers_Right,:);

%Calculating Coefficient of Pressure Along Right Symmetry:
c_p_Right = 1-(Right_mags(:,1).^2 + Right_mags(:,2).^2);

%Symmetry Left:
Left = Mesh.BC(6);
Left = cell2mat(Left);

%Computing Magnitude of Velocity along Right Symmetry:
Element_Numbers_Left = Left(:,3);
Left_mags = Gradient_w(Element_Numbers_Left,:);

%Calculating Coefficient of Pressure Along Right Symmetry:
c_p_Left = 1-(Left_mags(:,1).^2 + Left_mags(:,2).^2);


%Cslculating lift for each portion, and then summing all-together to
%determine lift coefficient:

%Determining Normals for Top:
n_Top = zeros(size(Top,1),2);
n_Top(:,2) = 1;
%Determining Normals for Right:
n_Right = zeros(size(Right,1),2);
n_Right(:,2) = 1;
%Determining Normals for Left:
n_Left = zeros(size(Left,1),2);
n_Left(:,2) = -1;


%Determining Sectional Lift for Top:
delta_p_Top = .5*1.225 *c_p_Top;
InnerX_Top = -delta_p_Top.*n_Top(:,1);
InnerY_Top = -delta_p_Top.*n_Top(:,2);

Inner_Top = [InnerX,InnerY];
for i=1:size(Inner_Top,1)
L_e_Top(i) = -dot(Inner_Top(i,:),Dists(i,:));
end

%Computing c_l for top boundary:
c_l_Top = trapz(L_e_Top)/((1/2)*1.225);

%Determining Sectional Lift for Right:
delta_p_Right = .5*1.225 *c_p_Right;
InnerX_Right = -delta_p_Right.*n_Right(:,1);
InnerY_Right = -delta_p_Right.*n_Right(:,2);

Inner_Right = [InnerX_Right,InnerY_Right];
for i=1:size(Inner_Right,1)
L_e_Right(i) = -dot(Inner_Right(i,:),Dists(i,:));
end

%Computing c_l for top boundary:
c_l_Right = trapz(L_e_Right)/((1/2)*1.225);

%Determining Sectional Lift for Left:
delta_p_Left = .5*1.225 *c_p_Left;
InnerX_Left = -delta_p_Left.*n_Left(:,1);
InnerY_Left = -delta_p_Left.*n_Left(:,2);

Inner_Left = [InnerX_Left,InnerY_Left];
for i=1:size(Inner_Left,1)
L_e_Left(i) = -dot(Inner_Left(i,:),Dists(i,:));
end

%Computing c_l for top boundary:
c_l_Left = trapz(L_e_Left)/((1/2)*1.225);


%Total c_l:
c_l_Total = c_l_Top + c_l_Right + c_l_Left



100 - 100*(c_l_Total/c_l)




%Found Values:
c_l = [0.5421 0.5458 0.5442 0.5401];
x = [868 3513 8655 36609]


%Plotting:
figure(6)
plot(x,c_l,'-b*')
title('Coefficient of Lift Based upon Number of Elements')
xlabel('Number of Elements')
ylabel('Coefficient of Lift')