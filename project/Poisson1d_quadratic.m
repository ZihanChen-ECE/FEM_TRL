%solving the 1D poisson equation with FEM.
%author: Zihan Chen
%modified time: March 23, 2016
%instruction: This program aims to solve the 1D poisson equation:
%                    -u''(x) = pi^2*sin(pi*x), x defines on [0,1]
%with Dirichlet BC as u(0) = u(1) = 0
%The analitical solution is: u(x) = sin(pi*x)
%This time I use quadratic element
%-----------------------------------------------------------------%

clear
clc

%%
%Space discretize
xl = 0;
xr = 1;
l = xr-xl;
N = 20;
h = l/N;
%number of elements. number of nodes will be N
u = zeros(N+1,1);
ue = zeros(3,1);
K = zeros(N+1,N+1);
Ke = zeros(3,3);
%find the discretize matrix in one element
alpha = 1;
beta = 0;
le = 2*h;
Ke(1,:) = [alpha*7/3/le+beta*2*le/15, -alpha*8/3/le+beta*le/15, alpha*1/3/le-beta*le/30];
Ke(2,:) = [-alpha*8/3/le+beta*le/15, alpha*16/3/le+beta*8*le/15, -alpha*8/3/le+beta*le/15];
Ke(3,:) = [alpha*1/3/le-beta*le/30, -alpha*8/3/le+beta*le/15, alpha*7/3/le+beta*2*le/15];
for i = 1:N-1
    K(i,i) = K(i,i) + Ke(1,1);
    K(i,i+1) = K(i,i+1) + Ke(1,2);
    K(i,i+2) = K(i,i+2) + Ke(1,3);
    K(i+1,i) = K(i+1,i) + Ke(2,1);
    K(i+1,i+1) = K(i+1,i+1) + Ke(2,2);
    K(i+1,i+2) = K(i+1,i+2) +Ke(2,3);
    K(i+2,i) = K(i+2,i) + Ke(3,1);
    K(i+2,i+1) = K(i+2,i+1) + Ke(3,2);
    K(i+2,i+2) = K(i+2,i+2) + Ke(3,3);
end

%apply the BC Dirichlet
Knew = K(2:N,2:N);
%find the source matrix    

%integrate be
he = h/N;
be = zeros(3,N-1);
xp = xl:h:xr;

f = @(x) pi^2*sin(pi*x);

be2 = zeros(3,N-1);

for j = 1:N-1 %element
xe1 = xp(j);
xe2 = xp(j+1);
xe3 = xp(j+2);

Ne1f = @(x) (x-xe2).*(x-xe3)/(xe1-xe2)/(xe1-xe3)*pi^2.*sin(pi*x);
Ne2f = @(x) (x-xe1).*(x-xe3)/(xe2-xe1)/(xe2-xe3)*pi^2.*sin(pi*x);
Ne3f = @(x) (x-xe1).*(x-xe2)/(xe3-xe1)/(xe3-xe2)*pi^2.*sin(pi*x);

be(1,j) = integral(Ne1f,xe1,xe3);
be(2,j) = integral(Ne2f,xe1,xe3);
be(3,j) = integral(Ne3f,xe1,xe3);

% fc =1/3*(f(xe1)+f(xe2)+f(xe3));
% be2(:,j) = [le/6*f(xe1);2/3*le*f(xe2);le/6*f(xe3)];
end

%norm (be(1,:)-be2(1,:))

% for j = 1:N-1
%     be(:,j) = [le/6;2/3*le;le/6];
% end
% b = zeros(N+1,1);
% b(1) = be(1,1);
% b(2) = be(2,1)+be(1,2);
% b(N) = be(2,N-1)+be(3,N-2);
% b(N+1) = be(3,N-1);
% for i = 3:N-1
%    b(i) = be(3,i-2)+be(2,i-1)+be(1,i);
% end
b = zeros(N+1,1);
for j = 1:N-1
   for k = 1:3
       b(j+k-1) = b(j+k-1) + be(k,j);
   end
end

bnew = b(2:end-1);

%%
%solving the matrix
u(2:end-1) = Knew\bnew;
U = sin(pi*xp);
plot(xp,u,'bo-')
hold on
plot(xp,U,'k--','linewidth',2)

%%
%finite difference 
B1 = ones(N-2,1);
B2 = ones(N-1,1);
B3 = B1;
A = diag(B1,-1)-2*diag(B2,0)+diag(B3,1);
A = -1/h^2*A;
bd = pi^2*sin(pi*xp(2:end-1));
bd = bd';
%bd = ones(N-1,1);
ud = zeros(N+1,1);
ud(2:end-1) = A\bd;
plot(xp,ud,'rO')
xlabel('x');
ylabel('u(x) or u_j');
title('1D FEM and FD for Poisson equation')
legend('U_F_E_M','U_t_r_u_e','U_F_D')
%error1 = norm(u-ud,1);
%error2 = norm(ud-U',1);

