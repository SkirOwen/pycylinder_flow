clc
clear all
% close all
%%
maxT   = 5000;  % total number of iterations
tPlot  = 1*50;      % cycles

% load('TestSymVars.mat')
load('TestSymVars.mat')

x = x-obst_x;
y = y-obst_y;
d = 2*obst_r;

%% visualize data
% close all
myContours = linspace(0, 1.4*uMax, 11);
for cycle = 2000:tPlot:maxT; %maxT
    %fileName = sprintf('Re100Data/velocity_%d.mat', cycle);
    fileName = sprintf('Testvelocity_%d.mat', cycle);
    load(fileName)
    u = sqrt(vx.^2 + vy.^2);

    figure(2), hold off
    contourf(x/d,y/d,u,myContours,'edgecolor','none'), hold on
    daspect([1 1 1])
    xlabel('x/d')
    ylabel('y/d')
    title(sprintf('Time step %d', cycle))
    %plot(x(300,50)/d,y(300,50)/d,'ko')
    colormap(viridisCMap)
    colorbar
    caxis([min(myContours),max(myContours)])
end

%% load data at a given point and plot
% k = 1;
% clear up
% for cycle = 5000:tPlot:maxT
%     fileName = sprintf('Re100Data/velocity_%d.mat', cycle);
%     load(fileName)
%     u = sqrt(vx.^2 + vy.^2);
%     up(k) = u(300,50);
%     k = k+1;
% end

%%
load('u_vs_t_downstream_v1.mat')
ti = (1:length(up));
figure(10), clf
plot(ti,up,'-r')
xlabel('Timestep')
ylabel('Flow speed, m/s')
grid on

