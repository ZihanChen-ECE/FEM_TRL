clear
clc


%%
%---------------------------------------------------------------------%
% Read and pre-condition the mesh file
%---------------------------------------------------------------------%
%gmsh_filename = 'example_2d.msh';
gmsh_filename = 'WG.msh';

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
%2D wave equation with FEM method
%---------------------------------------------------------------------%
 dt = 4e-10;
 omega = 2*pi*1.5*1e9;
 c = 3e8;
 kw = omega/c;
  t = 0:dt:4e-9;
 u = zeros(node_num,length(t));
 step = 1;
 for t = 0:dt:8e-9
    [Nelem,Nvertex] = size(element_nodeT);
    ENode = zeros(1,Nvertex);
    Ex = zeros(Nvertex,1);
    Ey = zeros(Nvertex,1);
    ae = zeros(Nvertex,1);
    be = zeros(Nvertex,1);
    ce = zeros(Nvertex,1);
    Eb = zeros(Nvertex,1);

    
    K = zeros(node_num,node_num);
    b = zeros(node_num,1);
  
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

        %calculate Kij for the element
        k(1,1) = (1/(4*EArea))*((Ey(2)-Ey(3))^2+(Ex(3)-Ex(2))^2)+-EArea/12*kw^2*(1+1);
        k(1,2) = (1/(4*EArea))*((Ey(2)-Ey(3))*(Ey(3)-Ey(1))+(Ex(3)-Ex(2))*(Ex(1)-Ex(3)))+-EArea/12*kw^2*(1);
        k(1,3) = (1/(4*EArea))*((Ey(2)-Ey(3))*(Ey(1)-Ey(2))+(Ex(3)-Ex(2))*(Ex(2)-Ex(1)))+-EArea/12*kw^2*(1);
        k(2,1) = (1/(4*EArea))*((Ey(3)-Ey(1))*(Ey(2)-Ey(3))+(Ex(1)-Ex(3))*(Ex(3)-Ex(2)))+-EArea/12*kw^2*(1);
        k(2,2) = (1/(4*EArea))*((Ey(3)-Ey(1))^2+(Ex(1)-Ex(3))^2)+-EArea/12*kw^2*(1+1);
        k(2,3) = (1/(4*EArea))*((Ey(3)-Ey(1))*(Ey(1)-Ey(2))+(Ex(1)-Ex(3))*(Ex(2)-Ex(1)))+-EArea/12*kw^2*(1);
        k(3,1) = (1/(4*EArea))*((Ey(1)-Ey(2))*(Ey(2)-Ey(3))+(Ex(2)-Ex(1))*(Ex(3)-Ex(2)))+-EArea/12*kw^2*(1);
        k(3,2) = (1/(4*EArea))*((Ey(1)-Ey(2))*(Ey(3)-Ey(1))+(Ex(2)-Ex(1))*(Ex(1)-Ex(3)))+-EArea/12*kw^2*(1);
        k(3,3) = (1/(4*EArea))*((Ey(1)-Ey(2))^2+(Ex(2)-Ex(1))^2)+-EArea/12*kw^2*(1+1);

        %update K with kij
        for j = 1:Nvertex
            for m = 1:Nvertex
            K(ENode(j),ENode(m)) = K(ENode(j),ENode(m)) + k(j,m);
            end
        end   

        %source term
 
        for j = 1:Nvertex
            if(Ex(j)==0)
                Eb(j) =-EArea/3*1;
                b(ENode(j)) = b(ENode(j)) + Eb(j);
            end 
        end
    end
  
 
 
 %Apply BC. Here I took Dirichelet BC at y = 0 and 0.15, while Neumann at x = 0 and 0.7
Bnode = zeros(1,1);
m = 1;

for i = 1:node_num
       if node_xT(i,2)==0 || node_xT(i,2)==0.15 %|| node_xT(i,2)==0 || node_xT(i,2)==1%Depends on the scope of the region
               Bnode(m) = i;
               m = m+1;
       end
end

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
u(:,step) = M\unew;
u(:,step) = u(:,step)*sin(omega*t);


step = step+1;
 end
 
for ss= 1:step-1 
figure;
trisurf(element_nodeT,node_xT(:,1),node_xT(:,2),u(:,ss));
colorbar; 
shading interp
axis([0 0.7 0 0.15])
xlabel('x');
ylabel('y');
title('2D FEM solution');
pause(0.5);
end
 
 
 
 
 
 
 
 
 