clc
clear
close all
tic()

%% preliminaries
ACC_Colors
m = 100;
n = m;
h = 2/m;
Int = linspace(-1,1,m+1);
XX = Int(1:(end-1));        % delta x = 0.02
u = 1;
U_init = [zeros(1,30),ones(1,10),zeros(1,60)];
time_init = 0;
delta_t = 0.00005;
num_steps = 250;
coords = [-1 1 -3 3];
r = 0;
C = u*delta_t/h;
%% 

%% Experiment

set(gcf, 'Position', [25, 25, 1900, 900])
plot(XX,U_init,'LineWidth',3,'Color',Color(:,11))
xlabel('x')
ylabel('u')
title('approximation of advection equation, C = 0.025')
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('Advection_1D.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

%% 4th order centered space, forward Euler in time

U = U_init;

for ii=1:num_steps
    U = U + C.*PDE_Finite_Diff_p(U,1,1,4,h);
    plot(XX,U','LineWidth',3,'Color',Color(:,11))
    xlabel('x')
    ylabel('u')
    title('approximation of advection equation, C = 0.025')
    axis(coords)
    drawnow()
    writeVideo(vidObj, getframe(gcf));
    
end
%%


close(vidObj);

toc()