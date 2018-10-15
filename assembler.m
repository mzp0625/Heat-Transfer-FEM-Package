function [K] = assembler(K_ele,conn,nn,ne)
K = zeros(nn,nn);
for i = 1:ne
    k11 = K_ele(1:2,1:2,i);
    k12 = K_ele(1:2,3:4,i);
    k21 = K_ele(3:4,1:2,i);
    k22 = K_ele(3:4,3:4,i);
    K(conn(i,1):conn(i,2),conn(i,1):conn(i,2)) =K(conn(i,1):conn(i,2),conn(i,1):conn(i,2))+k11;
    K(conn(i,1):conn(i,2),conn(i,3):conn(i,4)) =K(conn(i,1):conn(i,2),conn(i,3):conn(i,4))+k12;
    K(conn(i,3):conn(i,4),conn(i,1):conn(i,2)) =K(conn(i,3):conn(i,4),conn(i,1):conn(i,2))+k21;
    K(conn(i,3):conn(i,4),conn(i,3):conn(i,4)) =K(conn(i,3):conn(i,4),conn(i,3):conn(i,4))+k22;
end
end