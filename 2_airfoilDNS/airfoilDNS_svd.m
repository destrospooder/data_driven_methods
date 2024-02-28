%{
ME493 - Methods of Data-Driven Control
Benjamin Aziel

really needs refactoring
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
r = 4;
X_approx = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';

t_star = linspace(1, nt, 6); % timestep of interest (max r). will pick these better based on temporal amplitudes later

figure;
for k = 1:length(t_star)
    subplot(2, 3, k);
    ux_approx = reshape(X_approx(1:nx*ny, t_star(k)), nx, ny);
    contourf(x, y, ux_approx', linspace(-0.4, 0.2, 40));
    hold on;
    title(t_star(k))
    caxis([-0.4, 0.2]);
    colorbar;
    plot(xa(:,:),ya(:,:),'r-');  % plot all airfoil locations
end
plot(xa(:,:),ya(:,:),'k-')  % plot all airfoil locations
sgtitle(['ux reconstruction - rank 4']);

figure;
for k = 1:length(t_star)
    subplot(2, 3, k);
    uy_approx = reshape(X_approx(nx*ny+1:end, t_star(k)), nx, ny);
    contourf(x, y, uy_approx', linspace(-0.15, 0.15, 40));
    hold on;
    title(t_star(k))
    caxis([-0.15, 0.15]);
    colorbar;
    plot(xa(:,:),ya(:,:),'r-');  % plot all airfoil locations
end
sgtitle(['uy reconstruction - rank 4']);

%%
% DMD time

tv = 0; % truncation value

params_file = fullfile('airfoilDNS_parameters.h5');
dt_field = h5read(params_file, '/dt_field');

[Phi, Lambda, Atilde, Amplitudes] = calc_dmd(X, tv);

Eigscts = log(Lambda) / dt_field;

figure;
stem(imag(Eigscts)/(2*pi), abs(Amplitudes.*(Lambda))/max(abs(Amplitudes.*(Lambda))),'-','linewidth',2);
xlim([-0.02 1.5]);
ylim([0 1.2]);
xlabel('Frequency'), ylabel('DMD mode amplitude (scaled)');

%%
% SINDy time!

temp_amps = [V(:, 1:6) * S(1:6, 1:6)]; % x
t_field; % t

dtemp_amps = diff(temp_amps);
dtemp_amps = vertcat(dtemp_amps, dtemp_amps(400, :)); % to fill it out

polyorder = 2;
n = 6; % dofs
Theta = poolData(temp_amps, n, polyorder);
m = size(Theta,2);

lambda = 0.025;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta, dtemp_amps, lambda, n);
poolDataLIST({'x1', 'x2', 'x3', 'x4', 'x5', 'x6'}, Xi, n, polyorder);

%%
syms aleph
plot(t_field - 50, SV(1)*V(:, 1))
xlim([0 40])
hold on
fplot(0.0347/2 * aleph^2)
