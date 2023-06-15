C = load('Re100Data_Take2/u_vs_t_downstream_idx300_idx50_v1.mat');
F = load('Re100Data_FineGrid_Take2/u_vs_t_downstream_idx600_idx100_v1.mat');

C.ti = [1:length(C.up)];
F.ti = [1:length(F.up)];
%%
addpath('M:\Shared drives\AOE Lab7 Tunnel Staff\MattSzoke\Codes')
%%
figure(1), clf, hold on
% plot(C.ti(1:N:end),C.up(1:N:end),'-b.','displayname','Coarse grid')
subplot(2,1,1)
N = 1;
plot(0.5*C.ti(1:N:end),C.up(1:N:end),'-k.','displayname','Coarse grid, dT = 1')
grid on
xlim([20000 25000])
legend
ylabel('U, m/s')
xlabel('Timestep')

subplot(2,1,2)
N = 50;
plot(0.5*C.ti(1:N:end),C.up(1:N:end),'-k.','displayname','Coarse grid, dT = 50')
grid on
xlim([20000 25000])
legend
ylabel('U, m/s')
xlabel('Timestep')

%%
[ C.F, C.UP ] = myFFTb( C.up, 10, 1);
[ F.F, F.UP ] = myFFTb( F.up, 20, 1);

figure(2), clf, hold on
loglog(C.F, C.UP,'-b'), hold on
loglog(F.F, F.UP,'-k'), hold on
 xlim([0 0.006]*10)
 
 %%
N = 50;
[ C.F_L, C.UP_L ] = myFFTb( C.up(1:N:end), 1, 1);
[ C.F, C.UP ]     = myFFTb( C.up, N, 1);
% [ F.F, F.UP ] = myFFTb( F.up(1:K:end), 1, 1);

figure(2), clf, hold on
loglog(C.F, C.UP,'-b','displayname','FFT using every timestep'), hold on
loglog(C.F_L, C.UP_L,'-r','displayname','FFT using every N-th timestep'), hold on
%  xlim([0 0.006]*10)
title('Original (coarse) grid')
xlabel('Nondim. frequency')
ylabel('FFT$(u|_{x_0})$','interpreter','latex')
xlim([0.04 0.042])
grid on
legend
