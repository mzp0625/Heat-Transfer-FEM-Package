function [sigma_phi] = phi_sum(vector)
phi_1 = @(z1,z2) 1/4*(1-z1).*(1-z2).*(-1<=z1<=1).*(-1<=z2<=1);
phi_2 = @(z1,z2) 1/4*(1+z1).*(1-z2).*(-1<=z1<=1).*(-1<=z2<=1);
phi_3 = @(z1,z2) 1/4*(1-z1).*(1+z2).*(-1<=z1<=1).*(-1<=z2<=1);
phi_4 = @(z1,z2) 1/4*(1+z1).*(1+z2).*(-1<=z1<=1).*(-1<=z2<=1);
sigma_phi = @(z1,z2) vector(1)*phi_1(z1,z2)+vector(2)*phi_2(z1,z2)+vector(3)*phi_3(z1,z2)+vector(4)*phi_4(z1,z2);
end