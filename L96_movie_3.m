function L96_movie_3(Array_SqEnKF,Array_4DVar,Array_EDA,Array_En4DVar,Array_True,Array_Obs,...
    color1,color2,color3,color4,color5,color6,jump,ObsTimes)

ACC_Colors

coords = [-1.5 1.5 -1.5 1.5 -9 11];
ndim = size(Array_True,1);
mdim = size(Array_Obs,1);

theta = linspace(0,2*pi,ndim+1);
theta = theta(1:end-1);

X = cos(theta);
Y = sin(theta);
Xobs = X(1:2:end);
Yobs = Y(1:2:end);

exp_time = size(Array_True,2);

dasz = 17;
obssz1 = 20;
obssz2 = 10;
truesz = 3.5;
counter = 1;

figure()
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot3(X,Y,Array_SqEnKF(:,1),'.','Color',Color(:,color1),'MarkerSize',dasz);
hold on
h2 = plot3(X,Y,Array_4DVar(:,1),'.','Color',Color(:,color2),'MarkerSize',dasz);
h3 = plot3(X,Y,Array_EDA(:,1),'.','Color',Color(:,color3),'MarkerSize',dasz);
h4 = plot3(X,Y,Array_En4DVar(:,1),'.','Color',Color(:,color4),'MarkerSize',dasz);
h5 = plot3(X,Y,Array_True(:,1),'Color',Color(:,color5),'LineWidth',truesz);
h6 = plot3(Xobs,Yobs,Array_Obs(:,1),'*','Color',Color(:,color6),'MarkerSize',1);
plot3(Xobs,Yobs,Array_Obs(:,1),'.','Color',Color(:,color6),'MarkerSize',1);
xlabel('x')
ylabel('y')
zlabel('z')
axis(coords)
legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
hold off
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('L96_movie_3.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

for ii=2:exp_time
    if ii==ObsTimes(counter)
        h1 = plot3(X,Y,Array_SqEnKF(:,ii),'.','Color',Color(:,color1),'MarkerSize',dasz);
        hold on
        h2 = plot3(X,Y,Array_4DVar(:,ii),'.','Color',Color(:,color2),'MarkerSize',dasz);
        h3 = plot3(X,Y,Array_EDA(:,ii),'.','Color',Color(:,color3),'MarkerSize',dasz);
        h4 = plot3(X,Y,Array_En4DVar(:,ii),'.','Color',Color(:,color4),'MarkerSize',dasz);
        h5 = plot3(X,Y,Array_True(:,ii),'Color',Color(:,color5),'LineWidth',truesz);
        h6 = plot3(Xobs,Yobs,Array_Obs(:,counter),'*','Color',Color(:,color6),'MarkerSize',obssz2);
        plot3(Xobs,Yobs,Array_Obs(:,counter),'.','Color',Color(:,color6),'MarkerSize',obssz1);
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis(coords)
        legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    elseif ii==ObsTimes(counter)+1
        h1 = plot3(X,Y,Array_SqEnKF(:,ii),'.','Color',Color(:,color1),'MarkerSize',dasz);
        hold on
        h2 = plot3(X,Y,Array_4DVar(:,ii),'.','Color',Color(:,color2),'MarkerSize',dasz);
        h3 = plot3(X,Y,Array_EDA(:,ii),'.','Color',Color(:,color3),'MarkerSize',dasz);
        h4 = plot3(X,Y,Array_En4DVar(:,ii),'.','Color',Color(:,color4),'MarkerSize',dasz);
        h5 = plot3(X,Y,Array_True(:,ii),'Color',Color(:,color5),'LineWidth',truesz);
        h6 = plot3(Xobs,Yobs,Array_Obs(:,counter),'*','Color',Color(:,color6),'MarkerSize',obssz2);
        plot3(Xobs,Yobs,Array_Obs(:,counter),'.','Color',Color(:,color6),'MarkerSize',obssz1);
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis(coords)
        legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gcf));
        counter = counter + 1;
    else
        h1 = plot3(X,Y,Array_SqEnKF(:,ii),'.','Color',Color(:,color1),'MarkerSize',dasz);
        hold on
        h2 = plot3(X,Y,Array_4DVar(:,ii),'.','Color',Color(:,color2),'MarkerSize',dasz);
        h3 = plot3(X,Y,Array_EDA(:,ii),'.','Color',Color(:,color3),'MarkerSize',dasz);
        h4 = plot3(X,Y,Array_En4DVar(:,ii),'.','Color',Color(:,color4),'MarkerSize',dasz);
        h5 = plot3(X,Y,Array_True(:,ii),'Color',Color(:,color5),'LineWidth',truesz);
        h6 = plot3(Xobs,Yobs,Array_Obs(:,counter),'*','Color',Color(:,color6),'MarkerSize',0.1);
        plot3(Xobs,Yobs,Array_Obs(:,counter),'.','Color',Color(:,color6),'MarkerSize',0.1);
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis(coords)
        legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    end
end

close(vidObj);
end

