function [divergence_angle] = theta_angle(intensity, blabla)
    
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

    index = (x' >= max(x')/2);
    p = polyfit(x(index)',y(index)',1);
    yfit = p(2)+x'.*p(1)*1;
    divergence_angle = atan(p(1));
    
    if strcmp(blabla,'plot')
        figure,plot(x,y);
        hold on 
        plot(x,yfit);
        hold off
        xlabel ('z (\mum)');
        ylabel ('x (\mum)');
        title ('Divergence of Gaussian Beam');
    end
   
%     disp(divergence_angle)
end
