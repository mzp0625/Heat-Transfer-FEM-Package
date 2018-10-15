function [ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y)
% ne
ne = ne_x*ne_y;
% nn
nn_x = ne_x+1;
nn_y = ne_y+1;
nn = nn_x *nn_y;
% conn
temp_conn = zeros(ne_x,4,ne_y);
for i = 1: ne_x
    for j = 1: ne_y
        temp_conn(i,:,j) = [i+(j-1)*(nn_x) i+1+(j-1)*(nn_x) i+j*(nn_x) i+1+j*(nn_x)];
    end
end
conn = zeros(ne,4);
for i = 1:ne_x:ne
    conn(i:(i+ne_x-1),:) = temp_conn(:,:,ceil(i/ne_x));
end
% Loc_nodes
temp_loc = zeros(nn_x,2,nn_y);
x_dist = L(1)/ne_x;
y_dist  = L(2)/ne_y;
for i = 1:nn_x
    for j = 1:nn_y
        temp_loc(i,:,j) = [(i-1)*x_dist (j-1)*y_dist];
    end
end
for i = 1:nn_x:nn
    Loc_nodes(i:(i+nn_x-1),:) = temp_loc(:,:,ceil(i/nn_x));
end
% edge_conn
edge_1 = zeros(ne_y,2);
edge_2 = zeros(ne_x,2);
edge_3 = zeros(ne_y,2);
edge_4 = zeros(ne_x,2);
for i = 1:ne_y
    edge_1(i,:) = [nn_x+(i-1)*nn_x nn_x+i*nn_x];
    edge_3(i,:) = [1+(i-1)*nn_x 1+i*nn_x];
end
for i = 1:ne_x
    edge_2(i,:) = [i+ne_y*nn_x i+ne_y*nn_x+1 ];
    edge_4(i,:) = [i i+1];
end
end