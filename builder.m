function [K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne)
% location of nodes in X, y axis
X1 = zeros(ne,4);
X2 = zeros(ne,4);
for i = 1:ne
    ele_conn = conn(i,:);
    X1(i,:) = [Loc_nodes(ele_conn(1),1) Loc_nodes(ele_conn(2),1) Loc_nodes(ele_conn(3),1) Loc_nodes(ele_conn(4),1)];
    X2(i,:) = [Loc_nodes(ele_conn(1),2) Loc_nodes(ele_conn(2),2) Loc_nodes(ele_conn(3),2) Loc_nodes(ele_conn(4),2)];
end

% construct K_ele
K_ele = zeros(4,4,ne);
for i=1:ne
    % define phi_hat ,dphi_hat functions
    phi_1 = @(z1,z2) 1/4*(1-z1).*(1-z2);
    phi_2 = @(z1,z2) 1/4*(1+z1).*(1-z2);
    phi_3 = @(z1,z2) 1/4*(1+z1).*(1+z2);
    phi_4 = @(z1,z2) 1/4*(1-z1).*(1+z2);

    % define transformations
    x1_z = @(z1,z2) X1(i,1)*phi_1(z1,z2)+X1(i,2)*phi_2(z1,z2)+X1(i,3)*phi_3(z1,z2)+X1(i,4)*phi_4(z1,z2);
    x2_z = @(z1,z2) X2(i,1)*phi_1(z1,z2)+X2(i,2)*phi_2(z1,z2)+X2(i,3)*phi_3(z1,z2)+X2(i,4)*phi_4(z1,z2);

    % test location of center of master element on actual domain
    X_loc = x1_z(0,0);
    Y_loc = x2_z(0,0);
    if  h1<=X_loc && X_loc<=h1+d1 && h2<=Y_loc && Y_loc<=h2+d2
        k = k1;
    else
        k = k2;
    end
    % fill K_ele
    a = dim_x;
    b = dim_y;
    K_ele(1,1,i) = a^2+b^2;
    K_ele(1,2,i) = a^2/2-b^2;
    K_ele(1,3,i) = -a^2+b^2/2;
    K_ele(1,4,i) = -a^2/2-b^2/2;
    K_ele(2,3,i) = (-a^2-b^2)/2;
    K_ele(2,4,i) = -a^2+b^2/2;
    K_ele(3,4,i) = a^2/2-b^2;
    
    K_ele(2,2,i) = K_ele(1,1,i);
    K_ele(3,3,i) = K_ele(1,1,i);
    K_ele(4,4,i) = K_ele(1,1,i);
    K_ele(2,1,i) = K_ele(1,2,i);
    K_ele(3,1,i) = K_ele(1,3,i);
    K_ele(4,1,i) = K_ele(1,4,i);
    K_ele(3,2,i) = K_ele(2,3,i);
    K_ele(4,2,i) = K_ele(2,4,i);
    K_ele(4,3,i) = K_ele(3,4,i);
    
    K_ele(:,:,i) = k/(3*a*b)*K_ele(:,:,i);

end
end