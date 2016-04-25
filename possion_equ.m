%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function  :   ����Ԫ��Ⲵ�ɷ���(��������)
%%          :   ���ڲ�ͬ���⣬���޸�Դ���Լ��߽�������
%%See also  :   http://www.matlabsky.com/thread-10988-1-1.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
xLeft=0;                        %��
xRight=1;                       %��
yBottom=0;                      %��
yTop=1;                         %��

nNode=500;                      %��������
BigNum=1000;                    %�˴���
BoundN=floor(sqrt(nNode/2));    %ÿ���߽��ϵĵ���,���ֻ��һ������ֵ.floor��������ȡԼ��
%%%%%%%%%%%% �ڵ����� %%%%%%%%%%%%%
xNode=[linspace(xLeft,xRight,BoundN),xLeft*ones(1,BoundN),xRight*ones(1,BoundN),...
      linspace(xLeft,xRight,BoundN),rand(1,nNode-4*BoundN)];
yNode=[yBottom*ones(1,BoundN),linspace(yBottom,yTop,BoundN),linspace(yBottom,yTop,BoundN),...
     yTop*ones(1,BoundN),rand(1,nNode-4*BoundN)];  %��������ȥ���ظ��ĵ�
UsexyNode=unique([xNode(:),yNode(:)],'rows');
xNode=UsexyNode(:,1)';          yNode=UsexyNode(:,2)';
nNode=length(xNode);
%%%%%%%%%%%% ���������� %%%%%%%%%%%%%
Tri=delaunay(xNode,yNode);  %����������
TriN=size(Tri,1);           %�����θ���

GaussPointPosEta=[2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3];     GaussPointWeigh=[1/3,1/3,1/3];  %�������ڸ�˹���Լ�Ȩϵ��
Gausspoint2PosEta=[-1,1]/sqrt(3);                           GaussPoint2Weigh=[1,1];         %�߽���(��)���ֵ�Ȩ
xGaussPoint=xNode(Tri)*GaussPointPosEta';   %Gauss��x����
yGaussPoint=yNode(Tri)*GaussPointPosEta';   %Gauss��y����
%%%%%%%%%%%% ���ɷ�����Դ�� %%%%%%%%%%%%%
fGaussPoint=-2*pi*pi*sin(pi*xGaussPoint).*cos(pi*yGaussPoint);      %�޸ĵĵط�
%%%%%%%%%%%%% ����նȾ��� �� �����غ����� %%%%%%%%%%%%%%
TotalK=0;           TotalF=0;   
for TriI=1:TriN 
    TempA=[ones(3,1),xNode(Tri(TriI,:))',yNode(Tri(TriI,:))'];
    Ele_s=det(TempA)/2;                 %��Ԫ���
    TempA_1=inv(TempA); 
    
    Ele_f=GaussPointPosEta* ( fGaussPoint(TriI,:).*GaussPointWeigh *Ele_s )';       %���㵥Ԫ�غ�����
        TotalF=TotalF-accumarray(Tri(TriI,:)',Ele_f,[nNode,1],@sum);                %�ϳ������غ�����   
    
    Ele_k=TempA_1(2:3,:)'*TempA_1(2:3,:)*Ele_s; %���㵥Ԫ�նȾ���
        TriIndex=[Tri(TriI,[1,1,1,2,2,2,3,3,3])',Tri(TriI,[1,2,3,1,2,3,1,2,3])'];   %����նȾ����еı��
        TotalK=TotalK+accumarray(TriIndex,Ele_k(:),[nNode,nNode],@sum);             %�ϳ�����նȾ���       
end
%%%%%%%%%%%%% �߽��������� %%%%%%%%%%%%%%
%��һ�������
  index=find(xNode==0);   %��߽�Ľڵ�           %�޸ĵĵط�
  uAtB1=0;                %�޸ĵĵط�
  TotalIndex=sub2ind([nNode,nNode],index,index);      %�߽�ڵ�Ľڵ���
  TotalK(TotalIndex)=TotalK(TotalIndex)*BigNum;       %�ܸ���ӦԪ�س˴���
  TotalF(index)=TotalK(TotalIndex).*uAtB1;
        
  index=find(yNode==0);   %�±߽�Ľڵ�       %�޸ĵĵط�
  uAtB1=sin(pi*xNode(index));                 %�޸ĵĵط�
  TotalIndex=sub2ind([nNode,nNode],index,index);      %�߽�ڵ�Ľڵ���
  TotalK(TotalIndex)=TotalK(TotalIndex)*BigNum;       %�ܸ���ӦԪ�س˴���
  TotalF(index)=TotalK(TotalIndex).*uAtB1;    
    
  index=find(yNode==1);   %�ϱ߽�Ľڵ�       %�޸ĵĵط�
  uAtB1=(-sin(pi*xNode(index)));              %�޸ĵĵط�
  TotalIndex=sub2ind([nNode,nNode],index,index);      %�߽�ڵ�Ľڵ���
  TotalK(TotalIndex)=TotalK(TotalIndex)*BigNum;       %�ܸ���ӦԪ�س˴���
  TotalF(index)=TotalK(TotalIndex).*uAtB1;
%�ڶ���߽�����
  index=find(xNode==1);                           %�ұ߽�Ľڵ�           %�޸ĵĵط�
  CheckTri=ismember(Tri,index);                   %������Щ�ڵ�ĵ�Ԫ
  Bound2EleNum=find(sum(CheckTri,2)==2);          %�����ڵ�������(���������ڵ�)
 for TriI=1:length(Bound2EleNum)
     EleI2=Bound2EleNum(TriI);                   %�ڶ���߽��������漰�ĵ�Ԫ
     TriNodeIndex=find(CheckTri(EleI2,:)==1);    %�ߵ����˵��ڵ�Ԫ�ڵı��
     TotalNodeIndex=Tri(EleI2,TriNodeIndex);     %�ߵ����˵���ܱ��     
     LenEle=dist(yNode(TotalNodeIndex(2))-yNode(TotalNodeIndex(1)),...
     xNode(TotalNodeIndex(2))-xNode(TotalNodeIndex(1)));   %��Ԫ����
     Fia_i_g=[1-Gausspoint2PosEta;1+Gausspoint2PosEta]/2;        
     xyPos=[xNode(TotalNodeIndex);yNode(TotalNodeIndex)];        
     gPos=xyPos*Fia_i_g;                                         %��˹������
     xgPos=gPos(1,:);        ygPos=gPos(2,:);                    %��˹�������
     dudn=-pi*cos(pi*ygPos);         %�޸ĵĵط�
        
     EleF_Fix=GaussPoint2Weigh*LenEle/2.*dudn*Fia_i_g';          %�߻���
     TotalF(TotalNodeIndex)=TotalF(TotalNodeIndex)+EleF_Fix(:);  %���������غ�����
  end 
%%%%%%%%%%%% �ڵ�ĺ���ֵ %%%%%%%%%%%%
UNode=full(TotalK\TotalF);       
%%%%%%%%%%%%   ���ڳ�ͼ   %%%%%%%%%%%
U=sin(pi*xNode).*cos(pi*yNode);     %���۽�
subplot(2,2,1);
    trisurf(Tri,xNode,yNode,UNode);     axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
    xlabel('x');        ylabel('y');    title('����Ԫ��')
subplot(2,2,2); 
    trisurf(Tri,xNode,yNode,U(:));  axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
    xlabel('x');        ylabel('y');    title('���۽�')
subplot(2,2,3);     
    triplot(Tri,xNode,yNode);       axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
    xlabel('x');        ylabel('y');    title('����ֲ�')
subplot(2,2,4);
    trisurf(Tri,xNode,yNode,abs(UNode-U(:)));       axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
xlabel('x');        ylabel('y');    title('���ֲ�')
