function [divergence_angle,xx,yy] = theta_angle(intensity,z_mesh,x_mesh, wavelength,blabla)
    
    intensi_im = zeros(size(intensity));

    for i=1:length(intensity)
        max_ = max(intensity(:,i));
        threshold = max_/exp(2);
        intensi_im(:,i) = imbinarize(intensity(:,i),threshold);
    end
    
    intensi_im = intensi_im(1:size(intensi_im,1)/2,:);
%     figure,imshow(intensi_im)
    intensi_im_ed = edge(intensi_im); 
    intensi_im_ed = flipud(intensi_im_ed);
    [y, x] = find(intensi_im_ed == 1);
    y = y * x_mesh;
    x = x * z_mesh;
    
%     index = (x' >= max(x')/10*9);
%     p = polyfit(x(index)',y(index)',1);

    p = polyfit([0,mean([x(end),x(end-1)])],[0,mean([y(end),y(end-1)])],1);
    yfit = p(2)+x'.*p(1)*1;
    
    p1 = p(1);
    divergence_angle = atan(p1);
    divergence_angle_mrad = atan(p1)*1000;
    wo = min(y);
    zr = pi*wo^2/wavelength;
    R = max(x)*(1+(zr/max(x))^2);
    
    thetas = linspace(0,divergence_angle,100);
    xx = linspace(max(x)*cos(divergence_angle),max(x),100);
    yy = sqrt(max(x)^2-xx.^2);
    
%     curvature = linspace(R*cos(divergence_angle),R,100);
%     RR = R*sin(thetas);
    
    
    if strcmp(blabla,'plot')
        figure,plot(x,y);
        hold on 
        plot(x,yfit);
        plot(xx,yy);
        hold off
        
        xlabel ('z (\mum)');
        ylabel ('x (\mum)');
%         disp([num2str(max(x)),' ', num2str(max(y))])
        title ('Divergence of Gaussian Beam');
        
        text (min(x)+max(x)*0.1,max(y)-max(y)*0.1,['Divergence Angle = ',num2str(divergence_angle*180/pi),...
            ' ° = ',num2str(divergence_angle_mrad), ' mrad '])
        text (min(x)+max(x)*0.1,max(y)-max(y)*0.2,['Curvature(R) at z = ',num2str(max(x)),'\mum = ',num2str(R),' \mum '])
        text (min(x)+max(x)*0.1,max(y)-max(y)*0.3,['Rayleigh range (Zr) = ',num2str(zr),' \mum '])
        text (min(x)+max(x)*0.1,max(y)-max(y)*0.4,['Waist (wo) = ',num2str(wo),' \mum '])

    end
   
%     disp(divergence_angle)
end
