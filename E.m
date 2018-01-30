step = 0.1;
[x,y] = meshgrid(-2:step:2,-2:step:2);

pos_x = 0.0;
pos_y = 0.0;

q_x = x - pos_x;
q_y = y - pos_y;

E_x = q_x./(sqrt(q_x.^2+q_y.^2).^3);
E_y = q_y./(sqrt(q_x.^2+q_y.^2).^3);
E_abs = sqrt(E_x.^2+E_y.^2);
contourf(x,y,E_abs);
streamslice(x,y,E_x,E_y)




% Q1_x = -1;
% Q1_y = 0;
% Q1 = 1;
% 
% 
% step=1/16;
% [x,y] = meshgrid(-2:step:2,-2:step:2);
% 
% % charge 1
% r1_x=x-Q1_x;
% r1_y=y-Q1_y;
% 
% r1_3=(r1_x.^2+r1_y.^2).^(3/2);
% 
% E1_x=Q1*r1_x./r1_3;
% E1_y=Q1*r1_y./r1_3;
% 
% 
% E_abs=sqrt((E1_x).^2+(E1_y).^2);
% 
% [C,h] = contourf(x,y,log10(E_abs),50);
% shading interp
% 
% set(h,'LineColor','none')
% colormap pink
% caxis([-0.5,1.0])
% 
% hold on
% s=streamslice(x,y,E1_x+E2_x, E1_y+E2_y,2);
% set(s,'color','r')
% axis equal;
% axis([-2,2,-2,2]);
% 
% title(num2str(k));













% %% Set up some function. 
% % Sine between -2*pi and 2*pi.
% x = (10*-pi:0.1:10*pi)'; % Note the transpose.
% y = sin(x);
% fid = figure;
% set(fid,'Position',[0,0,1900,1600])
% hold on
% % The final plot.
% plot(x,y, '*');
%  
% 
% %% Set up the movie.
% writerObj = VideoWriter('out.avi'); % Name it.
% writerObj.FrameRate = 60; % How many frames per second.
% open(writerObj); 
%  
% for i=1:size(y)-500      
%     % We just use pause but pretend you have some really complicated thing here...
%     pause(0.1);
%     figure(fid); % Makes sure you use your desired frame.
%     plot(x(i),y(i),'or');
%  
%     %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
%         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
%     %end
%  
% end
% hold off
% close(writerObj); % Saves the movie.
