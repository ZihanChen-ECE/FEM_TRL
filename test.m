xLeft=0;
xRight=1;
PossionF=@(x)-pi^2*sin(pi*x);
%PossionF=@(x)x-x+1;


nNode=101;
BigNum=10000;
gausspos=[-1,1]./sqrt(3);
fia_gauss=[(1-gausspos)/2;(1+gausspos)/2];

K=zeros(nNode,nNode);
F=zeros(nNode,1);
x=linspace(xLeft,xRight,nNode);
dx=diff(x);

for ele=1:(nNode-1)
        EleLen=dx(ele);                %�߶γ���
        EleLen_1=1/EleLen;
        
        K(ele  ,ele  )=K(ele  ,ele  )+EleLen_1;
        K(ele  ,ele+1)=K(ele  ,ele+1)-EleLen_1;
        K(ele+1,ele  )=K(ele+1,ele  )-EleLen_1;
        K(ele+1,ele+1)=K(ele+1,ele+1)+EleLen_1;
        
        x_gauss  =x(ele)+EleLen*(1+gausspos)/2;
        f_gauss  =PossionF(x_gauss');
        term=fia_gauss*f_gauss*EleLen/2;
        
        F([ele,ele+1])=F([ele,ele+1])-term;
end


%%%%%%%%%%%%%%%%%%%
%%%%%%�߽�����%%%%%
%��һ��߽�����:��һ���ڵ�
F(1)=F(1)-pi;                        %��߽��ü�

%�ڶ���߽�����:���һ���ڵ�
K(nNode,nNode)=K(nNode,nNode)*BigNum;
F(nNode)=K(nNode,nNode)*0;

u=K\F;
plot(x,u,'r.',x,sin(pi*x))

legend('����ֵ','����ֵ')