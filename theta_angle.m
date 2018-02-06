function theta = theta_angle(intensity, z,x)

    intensi_im = zeros(size(intensity));
    for i=1:length(intensity)
        max_ = max(intensity(:,i));
        threshold = max_/exp(2);
        intensi_im(:,i) = imbinarize(intensity(:,i),threshold);
    end
        
    theta = intensi_im;
end
