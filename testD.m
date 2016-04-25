clear
clc

%%
%rand('state',0);
[x,y] = meshgrid(1:6,1:5);
tri = delaunay(x,y);
trisurf(tri,x,y,zeros(size(x)))
sum = 0;
for i = 1:length(tri(:,1))
    for j = 1:3
        if tri(i,j) == 17
            sum = sum+1;
        end
    end
end