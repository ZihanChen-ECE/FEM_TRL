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
xl = 0;
xr = 1;

l = xr-xl;
N = 18;
for in = 1:8
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
xp = xl:h:xr;
f = @(x) pi^2*sin(pi*x);
be = zeros(2,N);
for j = 1:N
    x1 = xp(j);
    x2 = xp(j+1);
    Ne1f = @(x) (x2-x)/h*pi^2.*sin(pi*x);
    Ne2f = @(x) (x-x1)/h*pi^2.*sin(pi*x);
    be(1,j) = integral(Ne1f,x1,x2);
    be(2,j) = integral(Ne2f,x1,x2);
end

%find the source matrix    
%  x = 0:h:1;
% b1e = @(x1,x2) pi*cos(pi*x1) - 1/h*(sin(pi*x2)-sin(pi*x1));
% b2e = @(x1,x2) -pi*cos(pi*x2) - 1/h*(sin(pi*x1)-sin(pi*x2));
% be = zeros(N,2);
% for i = 1:N   %go through all the elements
%     x1 = x(i);
%     x2 = x(i+1);
%     be(i,1) = b1e(x1,x2);
%     be(i,2) = b2e(x1,x2);
% end




b = zeros(N+1,1);
b(1,1) = be(1,1);
b(end,1) = be(2,end);
for i = 2:N
   b(i) = be(2,i-1)+be(1,i);
end
bnew = b(2:end-1);
%%
%solving the matrix
u(2:end-1) = Knew\bnew;
U = sin(pi*xp);
% plot(xp,u,'b*-')
% hold on
% plot(xp,U,'k--','linewidth',2)

%%
%finite difference 
B1 = ones(N-2,1);
B2 = ones(N-1,1);
B3 = B1;
A = diag(B1,-1)-2*diag(B2,0)+diag(B3,1);
A = -1/h^2*A;
bd = pi^2*sin(pi*xp(2:end-1));
bd = bd';
ud = zeros(N+1,1);
ud(2:end-1) = A\bd;
% plot(xp,ud,'rO')
% xlabel('x');
% ylabel('u(x) or u_j');
% title('1D FEM (linear element) and FD for Poisson equation')
% legend('U_F_E_M','U_t_r_u_e','U_F_D')
error1(in) = norm(ud-U',inf);
error2(in) = norm(ud-U',2);
%%
%Error analysis
he = h/N;
uref = @(x) sin(pi*x);
error = zeros(l/he+1,1);
ui = zeros(l/he+1,1);
xe = xl:he:xr;
for j = 1:N
    x1 = xp(j);
    x2 = xp(j+1);
    Ne1 = @(x) (x2-x)/h;
    Ne2 = @(x) (x-x1)/h;
    
    phi = @(x) (x2-x)/h*u(j)+(x-x1)/h*u(j+1);
    xs = x1:he:x2;
    for k = 1:length(xs)-1
        ui((j-1)*N+k) = phi(xs(k));
    end
    error((j-1)*N+1:j*N) = abs(phi(xe((j-1)*N+1:j*N))-uref(xe((j-1)*N+1:j*N))); 
    
end
% 
% figure(2);
% plot(xe, ui,'b-o');
% hold on
% plot(xe, uref(xe), 'r-*');
% xlabel('x(m)');
% ylabel('u');
% title('Linear interpolation of 1D fem solution vs Analytical solution');
% legend('1D fem','analytical solution');
% 
% figure(3);
% plot(xe,error,'linewidth',1);
% xlabel('x(m)');
% ylabel('error');
% title('Error of Linear interpolation');

L1error(in) = norm(error,inf);
L2error(in) = norm(error,2);
N = N*2;
end

for in = 1:7
   FDratio(in) = error1(in)/error1(in+1); 
end

for in = 1:7
   FDratio2(in) = error2(in)/error2(in+1); 
end

for in = 1:7
   LinearRatio(in) = L1error(in)/L1error(in+1); 
end

for in = 1:7
   LinearRatio2(in) = L2error(in)/L2error(in+1); 
end


