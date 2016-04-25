%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Function  :   有限元求解泊松方程(矩阵区域)
%%          :   对于不同问题，请修改源项以及边界条件项
%%See also  :   http://www.matlabsky.com/thread-10988-1-1.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
xLeft=0;                        %左
xRight=1;                       %右
yBottom=0;                      %下
yTop=1;                         %上

nNode=500;                      %初步估计
BigNum=1000;                    %乘大数
BoundN=floor(sqrt(nNode/2));    %每个边界上的点数,这个只是一个粗略值.floor函数用来取约数
%%%%%%%%%%%% 节点生成 %%%%%%%%%%%%%
xNode=[linspace(xLeft,xRight,BoundN),xLeft*ones(1,BoundN),xRight*ones(1,BoundN),...
      linspace(xLeft,xRight,BoundN),rand(1,nNode-4*BoundN)];
yNode=[yBottom*ones(1,BoundN),linspace(yBottom,yTop,BoundN),linspace(yBottom,yTop,BoundN),...
     yTop*ones(1,BoundN),rand(1,nNode-4*BoundN)];  %下面三行去掉重复的点
UsexyNode=unique([xNode(:),yNode(:)],'rows');
xNode=UsexyNode(:,1)';          yNode=UsexyNode(:,2)';
nNode=length(xNode);
%%%%%%%%%%%% 三角形生成 %%%%%%%%%%%%%
Tri=delaunay(xNode,yNode);  %三角形生成
TriN=size(Tri,1);           %三角形个数

GaussPointPosEta=[2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3];     GaussPointWeigh=[1/3,1/3,1/3];  %三角形内高斯点以及权系数
Gausspoint2PosEta=[-1,1]/sqrt(3);                           GaussPoint2Weigh=[1,1];         %边界上(线)积分的权
xGaussPoint=xNode(Tri)*GaussPointPosEta';   %Gauss点x坐标
yGaussPoint=yNode(Tri)*GaussPointPosEta';   %Gauss点y坐标
%%%%%%%%%%%% 泊松方程中源项 %%%%%%%%%%%%%
fGaussPoint=-2*pi*pi*sin(pi*xGaussPoint).*cos(pi*yGaussPoint);      %修改的地方
%%%%%%%%%%%%% 总体刚度矩阵 和 总体载荷向量 %%%%%%%%%%%%%%
TotalK=0;           TotalF=0;   
for TriI=1:TriN 
    TempA=[ones(3,1),xNode(Tri(TriI,:))',yNode(Tri(TriI,:))'];
    Ele_s=det(TempA)/2;                 %单元面积
    TempA_1=inv(TempA); 
    
    Ele_f=GaussPointPosEta* ( fGaussPoint(TriI,:).*GaussPointWeigh *Ele_s )';       %计算单元载荷向量
        TotalF=TotalF-accumarray(Tri(TriI,:)',Ele_f,[nNode,1],@sum);                %合成总体载荷向量   
    
    Ele_k=TempA_1(2:3,:)'*TempA_1(2:3,:)*Ele_s; %计算单元刚度矩阵
        TriIndex=[Tri(TriI,[1,1,1,2,2,2,3,3,3])',Tri(TriI,[1,2,3,1,2,3,1,2,3])'];   %总体刚度矩阵中的编号
        TotalK=TotalK+accumarray(TriIndex,Ele_k(:),[nNode,nNode],@sum);             %合成总体刚度矩阵       
end
%%%%%%%%%%%%% 边界条件处理 %%%%%%%%%%%%%%
%第一类边条件
  index=find(xNode==0);   %左边界的节点           %修改的地方
  uAtB1=0;                %修改的地方
  TotalIndex=sub2ind([nNode,nNode],index,index);      %边界节点的节点编号
  TotalK(TotalIndex)=TotalK(TotalIndex)*BigNum;       %总刚相应元素乘大数
  TotalF(index)=TotalK(TotalIndex).*uAtB1;
        
  index=find(yNode==0);   %下边界的节点       %修改的地方
  uAtB1=sin(pi*xNode(index));                 %修改的地方
  TotalIndex=sub2ind([nNode,nNode],index,index);      %边界节点的节点编号
  TotalK(TotalIndex)=TotalK(TotalIndex)*BigNum;       %总刚相应元素乘大数
  TotalF(index)=TotalK(TotalIndex).*uAtB1;    
    
  index=find(yNode==1);   %上边界的节点       %修改的地方
  uAtB1=(-sin(pi*xNode(index)));              %修改的地方
  TotalIndex=sub2ind([nNode,nNode],index,index);      %边界节点的节点编号
  TotalK(TotalIndex)=TotalK(TotalIndex)*BigNum;       %总刚相应元素乘大数
  TotalF(index)=TotalK(TotalIndex).*uAtB1;
%第二类边界条件
  index=find(xNode==1);                           %右边界的节点           %修改的地方
  CheckTri=ismember(Tri,index);                   %包含这些节点的单元
  Bound2EleNum=find(sum(CheckTri,2)==2);          %边所在的三角形(含有两个节点)
 for TriI=1:length(Bound2EleNum)
     EleI2=Bound2EleNum(TriI);                   %第二类边界条件所涉及的单元
     TriNodeIndex=find(CheckTri(EleI2,:)==1);    %边的两端点在单元内的编号
     TotalNodeIndex=Tri(EleI2,TriNodeIndex);     %边的两端点的总编号     
     LenEle=dist(yNode(TotalNodeIndex(2))-yNode(TotalNodeIndex(1)),...
     xNode(TotalNodeIndex(2))-xNode(TotalNodeIndex(1)));   %单元长度
     Fia_i_g=[1-Gausspoint2PosEta;1+Gausspoint2PosEta]/2;        
     xyPos=[xNode(TotalNodeIndex);yNode(TotalNodeIndex)];        
     gPos=xyPos*Fia_i_g;                                         %高斯点坐标
     xgPos=gPos(1,:);        ygPos=gPos(2,:);                    %高斯点的坐标
     dudn=-pi*cos(pi*ygPos);         %修改的地方
        
     EleF_Fix=GaussPoint2Weigh*LenEle/2.*dudn*Fia_i_g';          %线积分
     TotalF(TotalNodeIndex)=TotalF(TotalNodeIndex)+EleF_Fix(:);  %修正总体载荷向量
  end 
%%%%%%%%%%%% 节点的函数值 %%%%%%%%%%%%
UNode=full(TotalK\TotalF);       
%%%%%%%%%%%%   后期成图   %%%%%%%%%%%
U=sin(pi*xNode).*cos(pi*yNode);     %理论解
subplot(2,2,1);
    trisurf(Tri,xNode,yNode,UNode);     axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
    xlabel('x');        ylabel('y');    title('有限元解')
subplot(2,2,2); 
    trisurf(Tri,xNode,yNode,U(:));  axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
    xlabel('x');        ylabel('y');    title('理论解')
subplot(2,2,3);     
    triplot(Tri,xNode,yNode);       axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
    xlabel('x');        ylabel('y');    title('网格分布')
subplot(2,2,4);
    trisurf(Tri,xNode,yNode,abs(UNode-U(:)));       axis equal;     axis([xLeft,xRight,yBottom,yTop]);  colorbar;
xlabel('x');        ylabel('y');    title('误差分布')
