clear
clc

%%
%---------------------------------------------------------------------%
% Read and pre-condition the mesh file
%---------------------------------------------------------------------%
%gmsh_filename = 'example_2d.msh';
gmsh_filename = 'microstrip.msh';

  fprintf ( 1, '\n' );
  fprintf ( 1, 'GMSH_IO_TEST02:\n' );
  fprintf ( 1, '  Read data from a file.\n' );
%
%  Get the data size.
%
  [ node_num, m, element_num, element_order ] = gmsh_size_read ( gmsh_filename );
%
%  Print the sizes.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Node data read from file "%s"\n', gmsh_filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of nodes = %d\n', node_num );
  fprintf ( 1, '  Spatial dimension = %d\n', m );
  fprintf ( 1, '  Number of elements = %d\n', element_num );
  fprintf ( 1, '  Element order = %d\n', element_order );
%
%  Get the data.
%
   [ node_x, element_node ] = gmsh_data_read ( gmsh_filename, m, node_num, element_order, element_num );

    node_xT = node_x';
    element_nodeT = element_node'; 
    
     %discard those useless points in element_nodeT
   for i = 1:element_num
      if(element_nodeT(i,1)==0)
          break;
      end
   end
   element_nodeT(i:end,:)=[]; 
    
    
   %%
%---------------------------------------------------------------------%
%2D static Electrical problem (poison's equation) with FEM method
%---------------------------------------------------------------------%
epsilonair = 1;
epsilonsub = 3;
[Nelem,Nvertex] = size(element_nodeT);
ENode = zeros(1,Nvertex);
Ex = zeros(Nvertex,1);
Ey = zeros(Nvertex,1);
ae = zeros(Nvertex,1);
be = zeros(Nvertex,1);
ce = zeros(Nvertex,1);
Eb = zeros(Nvertex,1);

u = zeros(node_num,1);
K = zeros(node_num,node_num);
b = zeros(node_num,1);

k = zeros(3,3);
for i = 1:Nelem
    ENode(1,:) = element_nodeT(i,:);  %Find the node number for each element
    %get the position of x and y for each element node
    for j = 1:Nvertex
        Ex(j) = node_xT(ENode(j),1);
        Ey(j) = node_xT(ENode(j),2);   
    end
    for j = 1:Nvertex
        ae(1) = Ex(2)*Ey(3)-Ey(2)*Ex(3);
        ae(2) = Ex(3)*Ey(1)-Ey(3)*Ex(1);
        ae(3) = Ex(1)*Ey(2)-Ey(1)*Ex(2);
        be(1) = Ey(2)-Ey(3);
        be(2) = Ey(3)-Ey(1);
        be(3) = Ey(1)-Ey(2);
        ce(1) = Ex(3)-Ex(2);
        ce(2) = Ex(1)-Ex(3);
        ce(3) = Ex(2)-Ex(1);
    end
    EArea = abs(1/2*(Ex(2)*Ey(3)-Ey(2)*Ex(3)-Ex(1)*(Ey(3)-Ey(2))+Ey(1)*(Ex(3)-Ex(2))));
    
    
    if min(Ey(:))>=0.06 %%upper space
        epsilon = epsilonair;
    else
        epsilon = epsilonsub; %%lower space
    end
    %calculate Kij for the element
    k(1,1) = epsilon*(1/(4*EArea))*((Ey(2)-Ey(3))^2+(Ex(3)-Ex(2))^2);
    k(1,2) = epsilon*(1/(4*EArea))*((Ey(2)-Ey(3))*(Ey(3)-Ey(1))+(Ex(3)-Ex(2))*(Ex(1)-Ex(3)));
    k(1,3) = epsilon*(1/(4*EArea))*((Ey(2)-Ey(3))*(Ey(1)-Ey(2))+(Ex(3)-Ex(2))*(Ex(2)-Ex(1)));
    k(2,1) = epsilon*(1/(4*EArea))*((Ey(3)-Ey(1))*(Ey(2)-Ey(3))+(Ex(1)-Ex(3))*(Ex(3)-Ex(2)));
    k(2,2) = epsilon*(1/(4*EArea))*((Ey(3)-Ey(1))^2+(Ex(1)-Ex(3))^2);
    k(2,3) = epsilon*(1/(4*EArea))*((Ey(3)-Ey(1))*(Ey(1)-Ey(2))+(Ex(1)-Ex(3))*(Ex(2)-Ex(1)));
    k(3,1) = epsilon*(1/(4*EArea))*((Ey(1)-Ey(2))*(Ey(2)-Ey(3))+(Ex(2)-Ex(1))*(Ex(3)-Ex(2)));
    k(3,2) = epsilon*(1/(4*EArea))*((Ey(1)-Ey(2))*(Ey(3)-Ey(1))+(Ex(2)-Ex(1))*(Ex(1)-Ex(3)));
    k(3,3) = epsilon*(1/(4*EArea))*((Ey(1)-Ey(2))^2+(Ex(2)-Ex(1))^2);

    %update K with kij
    for j = 1:Nvertex
        for m = 1:Nvertex
        K(ENode(j),ENode(m)) = K(ENode(j),ENode(m)) + k(j,m);
        end
    end   
    
    %source term
%     for j = 1:Nvertex
%         Eb(j) = EArea/3*2*pi^2*sin(pi*Ex(j))*cos(pi*Ey(j));
%         b(ENode(j)) = b(ENode(j)) + Eb(j);
%     end   
end

%%
%Apply BC.
%First I deal with Dirichlet BC for phi = 1, as the source
%Find the position of these nodes
Snode = zeros(1,1);
Sm = 1;

%%Then I deal with Dirichlet BC for phi = 0
%Find the position of these nodes
Bnode = zeros(1,1);
Bm = 1;

for i = 1:node_num
       if node_xT(i,2)==0.3 || node_xT(i,1)==0.55 || node_xT(i,2)==0 %|| node_xT(i,2)==0 || node_xT(i,2)==1%Depends on the scope of the region
               Bnode(Bm) = i;
               Bm = Bm+1;
       end
       if (node_xT(i,2)==0.06 && node_xT(i,1)>=0.275 && node_xT(i,1)<0.33) ...
               ||(node_xT(i,2)==0.1 && node_xT(i,1)>=0.275 && node_xT(i,1)<0.33)...
               ||(node_xT(i,1)==0.33 && node_xT(i,2)>=0.06 && node_xT(i,2)<=0.1)
               Snode(Sm) = i;
               Sm = Sm+1;
       end
end

%Modify K and b with source node
for i=1:length(Snode)
    K(Snode(i),:)=0;
    K(Snode(i),Snode(i))=1;
    b(Snode(i)) = 1;
end

%reduce K and b by Dirichlet BC
%mapping matrix
M = eye(node_num,node_num);

 Z = length(Bnode);
% 
while Z > 0
    K(Bnode(Z),:)=[];
    K(:,Bnode(Z))=[];
    b(Bnode(Z)) = [];
    M(Bnode(Z),:) = [];
    Z = Z-1;
end
 unew = zeros(length(u)-length(Bnode),1);

 %solver u by K and b
 unew = K\b;
u = M\unew;

%plot z
figure(1)
trisurf(element_nodeT,node_xT(:,1),node_xT(:,2),u);
colorbar; 
axis([0.275 0.55 0 0.3])
shading interp
xlabel('x');
ylabel('y');
title('2D FEM solution');




