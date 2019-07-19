function [ Bzrp ] = B_from_input_magn_bkg_3d( bkg )
% Sums up psi and 3d-disturbed part of magnetic field
% bkg = read_ascot('input.magn_bkg');

% Reading the result:
% phivec = linspace(0.0, 2.0 * pi, bkg.nphi_per_sector)';
% B_r(R,phi,z) = Bzpr(z,R,phi,1) = interp3(bkg.R, bkg.z, phivec, Bzpr(:,:,:,1), R, z, phi);
% B_phi(R,phi,z) = Bzpr(z,R,phi,2) = interp3(bkg.R, bkg.z, phivec, Bzpr(:,:,:,2), R, z, phi);
% B_z(R,phi,z) = Bzpr(z,R,phi,3) = interp3(bkg.R, bkg.z, phivec, Bzpr(:,:,:,3), R, z, phi);

% plot field strength
% B = sqrt(Bzrp(:,:,:,1).^2 + Bzrp(:,:,:,2).^2 + Bzrp(:,:,:,3).^2);
% pcolor(bkg.R, bkg.z, B(:,:,1));

% plot components in Rz-plane, phi = 0. i = 1,2,3
% pcolor(bkg.R, bkg.z, Bzrp(:,:,1,i));

tol = 1e-5;

% psi part
bkg.psi = bkg.psi / (2.0 * pi);
sp = fnxtr(spaps({bkg.R, bkg.z}, bkg.psi', tol));
sp_dR = fnder(sp, [1 0]);% dpsi / dR
sp_dz = fnder(sp, [0 1]);% dpsi / dz

% Repeat R values in z direction
Rmat = repmat(bkg.R, 1, length(bkg.z));

BR = (-fnval(sp_dz, {bkg.R, bkg.z}) ./ Rmat)';
Bz = (fnval(sp_dR, {bkg.R, bkg.z}) ./ Rmat)';
Bphi = bkg.Bphi;

BR = repmat(BR, [1, 1, bkg.nphi_per_sector]);
Bz = repmat(Bz, [1, 1, bkg.nphi_per_sector]);

% Add 3d part
BR = BR + bkg.BR;
Bz = Bz + bkg.Bz;   

Bzrp(:,:,:,1) = BR;
Bzrp(:,:,:,2) = Bphi;
Bzrp(:,:,:,3) = Bz;
end
