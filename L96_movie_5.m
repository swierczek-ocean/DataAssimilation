function L96_movie_5(Array1,Array2,color1,color2)

ACC_Colors

coords = [-1.5 1.5 -1.5 1.5 -9 11];
ndim = size(Array1,1);

theta = linspace(0,2*pi,ndim+1);
theta = theta(1:end-1);

X = cos(theta);
Y = sin(theta);

exp_time = size(Array1,2);

truesz1 = 5.5;
truesz2 = 3.2;

figure()
set(gcf, 'Position', [25, 25, 1600, 900])
plot3(X,Y,Array1(:,1),'.','Color',Color(:,color1),'MarkerSize',truesz1);
hold on
plot3(X,Y,Array2(:,1),'.','Color',Color(:,color2),'MarkerSize',truesz2);
xlabel('x')
ylabel('y')
zlabel('z')
axis(coords)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('L96_movie_5.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 45;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

for ii=2:exp_time
    plot3(X,Y,Array1(:,ii),'Color',Color(:,color1),'LineWidth',truesz1);
    hold on
    plot3(X,Y,Array2(:,ii),'Color',Color(:,color2),'LineWidth',truesz2);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis(coords)
    hold off
    drawnow()
    writeVideo(vidObj, getframe(gcf));
end

close(vidObj);
end

