%solving the 1D poisson equation with FEM.
%author: Zihan Chen
%modified time: March 23, 2016
%instruction: This program aims to solve the 1D poisson equation:
%                    -u''(x) = pi^2*sin(pi*x), x defines on [0,1]
%with Dirichlet BC as u(0) = u(1) = 0
%The analitical solution is: u(x) = sin(pi*x)
%-------------------------------------------------
clear
clc

%%
%Space discretize
l = 1;
N = 20;
h = l/N;
%number of elements. number of nodes will be N
u = zeros(N+1,1);
ue = zeros(2,1);
K = zeros(N+1,N+1);
Ke = zeros(2,2);
%find the discretize matrix in one element
Ke = [1/h,-1/h;-1/h,1/h];
for i = 1:N
    K(i,i) = K(i,i) + Ke(1,1);
    K(i,i+1) = K(i,i+1) + Ke(1,2);
    K(i+1,i) = K(i+1,i) + Ke(2,1);
    K(i+1,i+1) = K(i+1,i+1) + Ke(2,2);
end
%apply the BC
Knew = K(2:N,2:N);
%find the source matrix    
 x = 0:h:1;
b1e = @(x1,x2) pi*cos(pi*x1) - 1/h*(sin(pi*x2)-sin(pi*x1));
b2e = @(x1,x2) -pi*cos(pi*x2) - 1/h*(sin(pi*x1)-sin(pi*x2));
be = zeros(N,2);
for i = 1:N   %go through all the elements
    x1 = x(i);
    x2 = x(i+1);
    be(i,1) = b1e(x1,x2);
    be(i,2) = b2e(x1,x2);
end
b = zeros(N+1,1);
b(1,1) = be(1,1);
b(end,1) = be(end,2);
for i = 2:N
   b(i) = be(i-1,2)+be(i,1);
end
bnew = b(2:end-1);
%%
%solving the matrix
u(2:end-1) = Knew\bnew;
U = sin(pi*x);
plot(x,u,'bo-')
hold on
plot(x,U,'k--','linewidth',2)

%%
%finite difference 
B1 = ones(N-2,1);
B2 = ones(N-1,1);
B3 = B1;
A = diag(B1,-1)-2*diag(B2,0)+diag(B3,1);
A = -1/h^2*A;
bd = pi^2*sin(pi*x(2:end-1));
bd = bd';
ud = zeros(N+1,1);
ud(2:end-1) = A\bd;
plot(x,ud,'rO')
xlabel('x');
ylabel('u(x) or u_j');
title('1D FEM and FD for Poisson equation')
legend('U_F_E_M','U_t_r_u_e','U_F_D')
error1 = norm(u-U',1);
error2 = norm(ud-U',1);



