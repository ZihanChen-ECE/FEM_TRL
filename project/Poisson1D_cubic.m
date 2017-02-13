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
Ke = zeros(4,4);
%find the discretize matrix in one element
alpha = 1;
beta = 0;
le = 3*h;
Ke(1,:) = [alpha*37/10/le+beta*8*le/105, -alpha*189/40/le+beta*33*le/560, alpha*27/20/le-beta*3*le/140, -alpha*13/40/le+beta*19*le/1680];
Ke(2,:) = [-alpha*189/40/le+beta*33*le/560, alpha*54/5/le+beta*27*le/70, -alpha*297/40/le-beta*le*27/560, alpha*27/20/le-beta*3*le/140];
Ke(3,:) = [alpha*27/20/le-beta*3*le/140, -alpha*297/40/le-beta*le*27/560, alpha*54/5/le+beta*27*le/70, -alpha*189/40/le+beta*33*le/560];
Ke(4,:) = [-alpha*13/40/le+beta*19*le/1680, alpha*27/20/le-beta*3*le/140, -alpha*189/40/le+beta*33*le/560, alpha*37/10/le+beta*8*le/105];
for i = 1:N-2
    K(i,i) = K(i,i) + Ke(1,1);
    K(i,i+1) = K(i,i+1) + Ke(1,2);
    K(i,i+2) = K(i,i+2) + Ke(1,3);
    K(i,i+3) = K(i,i+3) + Ke(1,4);
    
    K(i+1,i) = K(i+1,i) + Ke(2,1);
    K(i+1,i+1) = K(i+1,i+1) + Ke(2,2);
    K(i+1,i+2) = K(i+1,i+2) + Ke(2,3);
    K(i+1,i+2) = K(i+1,i+3) + Ke(2,4);
    
    K(i+2,i) = K(i+2,i) + Ke(3,1);
    K(i+2,i+1) = K(i+2,i+1) + Ke(3,2);
    K(i+2,i+2) = K(i+2,i+2) + Ke(3,3);
    K(i+2,i+3) = K(i+2,i+3) + Ke(3,4);
    
    K(i+3,i) = K(i+3,i) + Ke(4,1);
    K(i+3,i+1) = K(i+3,i+1) + Ke(4,2);
    K(i+3,i+2) = K(i+3,i+2) + Ke(4,3);
    K(i+3,i+3) = K(i+3,i+3) + Ke(4,4);
    
end

%apply the BC Dirichlet
Knew = K(2:N,2:N);
%find the source matrix    

%integrate be
he = h/N;
be = zeros(4,N-2);
xp = xl:h:xr;

f = @(x) pi^2*sin(pi*x);

% be2 = zeros(4,N-2);

for j = 1:N-2 %element
xe1 = xp(j);
xe2 = xp(j+1);
xe3 = xp(j+2);
xe4 = xp(j+3);

Ne1f = @(x) (x-xe2).*(x-xe3).*(x-xe4)/(xe1-xe2)/(xe1-xe3)/(xe1-xe4)*pi^2.*sin(pi*x);
Ne2f = @(x) (x-xe1).*(x-xe3).*(x-xe4)/(xe2-xe1)/(xe2-xe3)/(xe2-xe4)*pi^2.*sin(pi*x);
Ne3f = @(x) (x-xe1).*(x-xe2).*(x-xe4)/(xe3-xe1)/(xe3-xe2)/(xe3-xe4)*pi^2.*sin(pi*x);
Ne4f = @(x) (x-xe1).*(x-xe2).*(x-xe3)/(xe4-xe1)/(xe4-xe2)/(xe4-xe3)*pi^2.*sin(pi*x); 

be(1,j) = integral(Ne1f,xe1,xe4);
be(2,j) = integral(Ne2f,xe1,xe4);
be(3,j) = integral(Ne3f,xe1,xe4);
be(4,j) = integral(Ne4f,xe1,xe4);

% fc =1/3*(f(xe1)+f(xe2)+f(xe3));
% be2(:,j) = [le/6*f(xe1);2/3*le*f(xe2);le/6*f(xe3)];
end

%norm (be(1,:)-be2(1,:))

% for j = 1:N-1
%     be(:,j) = [le/6;2/3*le;le/6];
% end
b = zeros(N+1,1);
for j = 1:N-2
   for k = 1:4
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

