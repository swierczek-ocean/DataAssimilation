function L96_movie_4(Array_1,Array_2,color1,color2,dim1,dim2,dim3,coords)

ACC_Colors

Array1 = Array_1([dim1,dim2,dim3],:);
Array2 = Array_2([dim1,dim2,dim3],:);

exp_time = size(Array1,2);

truesz1 = 5.5;
truesz2 = 3.2;

figure()
set(gcf, 'Position', [25, 25, 1600, 900])
plot3(Array1(1,1),Array1(2,1),Array1(3,1),'.','Color',Color(:,color1),'MarkerSize',truesz1);
hold on
plot3(Array2(1,1),Array2(2,1),Array2(3,1),'.','Color',Color(:,color2),'MarkerSize',truesz2);
xlabel('x')
ylabel('y')
zlabel('z')
axis(coords)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('L96_movie_4.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 45;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

for ii=2:floor(exp_time/2)-10
    plot3(Array1(1,1:2*ii),Array1(2,1:2*ii),Array1(3,1:2*ii),'Color',Color(:,color1),'LineWidth',truesz1);
    hold on
    plot3(Array2(1,1:2*ii),Array2(2,1:2*ii),Array2(3,1:2*ii),'Color',Color(:,color2),'LineWidth',truesz2);
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

