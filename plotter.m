function [] = plotter(a_j,ne_x,ne_y,Loc_nodes)
    X = Loc_nodes(1:ne_x+1,1);
    Y = Loc_nodes(1:(ne_x+1):length(Loc_nodes),2);
    a = zeros(ne_x+1,ne_y+1);
    nn_y = ne_y+1;
    for j = 1:nn_y
        a(:,j) = a_j(j+(j-1)*ne_x:j+(j)*ne_x);
    end
    a = a';
    figure;
    surf(X,Y,a);
    colorbar;
    xlabel('X coordinate');
    ylabel('Y coordinate');
    zlabel('T');
end
