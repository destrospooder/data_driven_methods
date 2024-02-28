%{

ME493 - Methods of Data-Driven Control
POD analysis of Low-Reynolds-number pitching airfoil 
    direct numerical simulations
Benjamin Aziel

%}

clear; clc;

grid_file = fullfile('airfoilDNS_grid.h5');
data_file = fullfile('airfoilDNS_a25f0p05.h5');

mean_correction = true;

x = h5read(grid_file, '/x');
y = h5read(grid_file, '/y');
nx = length(x);
ny = length(y);

t_field = h5read(data_file, '/t_field');
t_force = h5read(data_file, '/t_force');
nt = length(t_field);

ux = h5read(data_file, '/ux');
ux_reshaped = reshape(ux, nx * ny, nt);
uy = h5read(data_file, '/uy');
uy_reshaped = reshape(uy, nx * ny, nt);

xa = h5read(data_file, '/xa');
ya = h5read(data_file, '/ya');

X = [ux_reshaped; uy_reshaped];

if mean_correction
    X_mean = mean(X, 2);
    X = X - X_mean * ones(1,nt);
end

[U, S, V] = svd(X, 'econ');

%% 
SV = diag(S);

figure;
semilogy(linspace(1, length(S), length(S)), SV.^2, linewidth = 2);
xlabel('Indices');
ylabel('Squared Singular Values');
title('Squared SVs - No Editing')

figure;
semilogy(linspace(1, length(S)-1, length(S)-1), SV(1:end-1).^2, linewidth = 2);
xlabel('Indices');
ylabel('Squared Singular Values');
title('Squared SVs - First 400')

U_ux = U(1:length(ux_reshaped), :);
U_uy = U(length(uy_reshaped)+1:end, :);

MM = 0.01;
v = -1:0.1:1;
v(11)=[];

%%
figure;
for k = 1:6
    subplot(2, 3, k);
    contourf(x, y, transpose(reshape(U_ux(:,k), nx, ny)), MM * v);
    caxis([-MM MM]);
    colorbar;
end
sgtitle('ux spatial modes');

figure;
for k = 1:6
    subplot(2, 3, k);
    contourf(x, y, transpose(reshape(U_uy(:,k), nx, ny)), MM * v);
    caxis([-MM MM]);
    colorbar;
end
sgtitle('uy spatial modes');


figure;
for k = 1:6
    subplot(2, 3, k);
    plot(t_field, SV(k)*V(:, k));
    hold on;
end
sgtitle(['temporal amplitudes']);

%%

% Let's attempt reconstruction for rank r
r = 125;
X_approx = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';

t_star = 200; % timestep of interest (max 401)
ux_approx = reshape(X_approx(1:nx*ny, t_star), nx, ny);

figure;
contourf(x, y, ux_approx', linspace(-0.4, 0.2, 120));
caxis;
colorbar;
