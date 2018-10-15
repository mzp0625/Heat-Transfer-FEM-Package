clear all
%% 3.1 Shape Functions
z1 = -1:0.05:1;
z2 = z1;
vector_1 = [1 0 2 0];
vector_2 = [1 4 2 0];
sigma_phi_1 = phi_sum(vector_1);
sigma_phi_2 = phi_sum(vector_2);
sigma_phi_1_val = zeros(length(z1),length(z2));
sigma_phi_2_val = zeros(length(z1),length(z2));
for i = 1:length(z1)
    for j = 1:length(z2)
        sigma_phi_1_val(i,j) = sigma_phi_1(z1(i),z2(j));
        sigma_phi_2_val(i,j) = sigma_phi_2(z1(i),z2(j));
    end
end

figure
surf(z1,z2,sigma_phi_1_val);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
colorbar;
xlabel('z1');
ylabel('z2');
zlabel('sigma(a_iphi_i)');
title('sum of phi functions, v = [1 0 2 0]');

figure
surf(z1,z2,sigma_phi_2_val);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
colorbar;
xlabel('z1');
ylabel('z2');
zlabel('sigma(a_iphi_i)');
title('sum of phi functions, v = [1 4 2 0]');

%% 3.3 Mesher

%% 3.3.1
ne_x = 2;
ne_y = 2;
L =[0.75 0.5];
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
disp('the number of elements is');
disp(ne);
disp('the number of nodes is');
disp(nn);
disp('the element connectivity matrix is');
disp(conn);
disp('the node location matrix is ');
disp(Loc_nodes);
disp('the edge connectivity matrices are');
disp(edge_1);
disp(edge_2);
disp(edge_3);
disp(edge_4);

%% 3.3.2
k1 = 2;
k2 = 2;
d1 = 1;
d2 = 1;
h1 = 0;
h2 = 0;
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
disp('the first elemental stiffness matrix is');
disp(K_ele(:,:,1));

%% 3.4 Assembler
% nothing here
[K] = assembler(K_ele,conn,nn,ne);
%% 3.5.1
phi = -20;
R = zeros(nn,1);
[R] = flux(R,Loc_nodes,edge_2,phi);
disp('the load vector due to flux =1 on boundary 2 is');
disp(R);

%% 3.5.2
T = 50;
[R,K] = fixed(R,K,edge_4,T);
disp('the load vector due to T = 30 on boundary 4 is');
disp(R);

%% 3.6 what about other edges?
% nothing 

%% 3.7 solver
[a_j] = solver(K,R);

%% 3.8 plotter
plotter(a_j,ne_x,ne_y,Loc_nodes);

%% 4 Checking code
clear all
% numerical, case 1
L = [3 5];
h = [1 2];
d = [1 1];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 5;
k2 = 5;
ne_x = 10;
ne_y = 10;
T = 32;
phi = 10;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R] = flux(R,Loc_nodes,edge_3,phi);
[R,K] = fixed(R,K,edge_1,T);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('numerical, case 1');
% analytical, case 1
f = @(x) phi*x/k1+T;
figure;
fplot(f,[0 3]);
xlabel('edge 3 -X coordinate - edge 1');
ylabel('T');
title('analytical, case 1');
% numerical, case 2
T = 50;
phi = -20;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R] = flux(R,Loc_nodes,edge_2,phi);
[R,K] = fixed(R,K,edge_4,T);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('numerical, case 2');
% analytical, case 2
f = @(x) phi*x/k1+T;
figure;
fplot(f,[0 5]);
xlabel('edge 4 - Y coordinate - edge 2');
ylabel('T');
title('analytical, case 2');
%% 5 Simulation
%% simulation 1
clear all
L = [1 2];
h = [0.3 0.5];
d = [0.4 1];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 100;
k2 = 5;
ne_x = 10;
ne_y = 20;
T1 = 32;
T3 = 50;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R,K] = fixed(R,K,edge_1,T1);
[R,K] = fixed(R,K,edge_3,T3);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('Simulation 1');
%% simulation 2
clear all
L = [1 2];
h = [0.3 0.5];
d = [0.4 1];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 5;
k2 = 100;
ne_x = 10;
ne_y = 20;
T1 = 32;
T3 = 50;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R,K] = fixed(R,K,edge_1,T1);
[R,K] = fixed(R,K,edge_3,T3);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('Simulation 2');
%% simulation 3
clear all
L = [1 2];
h = [0.3 0.5];
d = [0.4 1];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 100;
k2 = 5;
ne_x = 10;
ne_y = 20;
T1 = 32;
T3 = 50;
phi2 = 100;
phi4 = -200;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R] = flux(R,Loc_nodes,edge_2,phi2);
[R] = flux(R,Loc_nodes,edge_4,phi4);
[R,K] = fixed(R,K,edge_1,T1);
[R,K] = fixed(R,K,edge_3,T3);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('Simulation 3');
%% simulation 4
clear all
L = [1 2];
h = [0.3 0.5];
d = [0.4 1];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 5;
k2 = 5;
ne_x = 10;
ne_y = 20;
T1 = 32;
phi2 = 100;
phi3 = 50;
phi4 = -200;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R] = flux(R,Loc_nodes,edge_2,phi2);
[R] = flux(R,Loc_nodes,edge_3,phi3);
[R] = flux(R,Loc_nodes,edge_4,phi4);
[R,K] = fixed(R,K,edge_1,T1);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('Simulation 4');
%% simulation 5
clear all
L = [1 2];
h = [0 0];
d = [0.6 1.4];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 5000;
k2 = 5;
ne_x = 10;
ne_y = 20;
T1 = 32;
T3 = 50;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R,K] = fixed(R,K,edge_1,T1);
[R,K] = fixed(R,K,edge_3,T3);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('Simulation 5');
%% simulation 6
clear all
L = [1 2];
h = [0 0];
d = [0.6 1.4];
h1 = h(1);
h2 = h(2);
d1 = d(1);
d2 = d(2);
k1 = 5;
k2 = 50;
ne_x = 20;
ne_y = 20;
T1 = 32;
T3 = 50;
phi4 = 500;
[ne,nn ,conn,Loc_nodes,edge_1,edge_2,edge_3,edge_4] = mesher(L,ne_x,ne_y);
dim_x = L(1)/ne_x;
dim_y = L(2)/ne_y;
[K_ele] = builder(dim_x,dim_y,k1,k2,h1,h2,d1,d2,Loc_nodes,conn,ne);
[K] = assembler(K_ele,conn,nn,ne);
R = zeros(nn,1);
[R] = flux(R,Loc_nodes,edge_4,phi4);
[R,K] = fixed(R,K,edge_1,T1);
[R,K] = fixed(R,K,edge_3,T3);
[a_j] = solver(K,R);
plotter(a_j,ne_x,ne_y,Loc_nodes);
title('Simulation 5');