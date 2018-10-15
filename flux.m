function[R] = flux(R,Loc_nodes,edge_array,phi)
w8 = [0.3626837833783620 0.3626837833783620 0.3137066458778873 0.3137066458778873 0.2223810344533745 0.2223810344533745 0.1012285362903763 0.1012285362903763];
z8 = [-0.1834346424956498 0.1834346424956498 -0.5255324099163290 0.5255324099163290 -0.7966664774136267 0.7966664774136267 -0.9602898564975363 0.9602898564975363];
j = 1:8;
% define phi in parametric space
phi_1 = @(z) (1-z)./2;
phi_2 = @(z) (1+z)./2;
dphi_1 = @(z) -1./2.*z.^0;
dphi_2 = @(z) 1./2.*z.^0;

edge_conn = [edge_array(:,1)' edge_array(end,end)];
if range(Loc_nodes(edge_conn,1))==0 % vertical line
    X = Loc_nodes(edge_conn,2);
elseif range(Loc_nodes(edge_conn,2)) ==0 % horizontal line
    X = Loc_nodes(edge_conn,1);
end

R_ele = zeros(2,1,length(X)-1);
for i = 1:length(X)-1
    J = @(z) X(i)*dphi_1(z)+X(i+1)*dphi_2(z);
    r_1 = @(z) phi*phi_1(z).*J(z);
    r_2 = @(z) phi*phi_2(z).*J(z);
    R_ele(1,1,i) = sum(w8(j).*r_1(z8(j)));
    R_ele(2,1,i) = sum(w8(j).*r_2(z8(j)));
end

for i = 1:length(X)-1
    R(edge_array(i,1))= R(edge_array(i,1))+R_ele(1);
    R(edge_array(i,2))= R(edge_array(i,2))+R_ele(2);
    
end
end