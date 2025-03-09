img = imread('C:\Users\giris\Downloads\keys.pnm');
sigma = 2;
low_threshold = 0.001;
high_threshold = 0.5;
subplot(2,2,1), imshow(img), title('Original Image');
edges = canny_edge_detection(img, sigma, low_threshold, high_threshold);
subplot(2,2,2), imshow(edges), title('Edge Map');
img_edges = img;
img_edges(edges==1) = 255;
subplot(2,2,3), imshow(img_edges), title('Edges Superimposed');
function edges = canny_edge_detection(img, sigma, ~, ~)
    filter_size = 2*ceil(3*sigma)+1;
    filter = fspecial('gaussian', filter_size, sigma);
    img_smoothed = conv2(img, filter, 'same');
    dx_filter = [-1 0 1; -sqrt(2) 0 sqrt(2); -1 0 1]/(2*sqrt(2));
    dy_filter = dx_filter';
    dx = conv2(img_smoothed, dx_filter, 'same');
    dy = conv2(img_smoothed, dy_filter, 'same');
    mag = sqrt(dx.^2 + dy.^2);
    theta = atan2d(dy, dx);
    [m, n] = size(img);
    edges = zeros(m, n);
    for i = 2:m-1
        for j = 2:n-1
            if (theta(i,j) < -22.5) || (theta(i,j) > 157.5)
                orientation = 0;
            elseif (theta(i,j) >= -22.5) && (theta(i,j) < 22.5)
                orientation = 1;
            elseif (theta(i,j) >= 22.5) && (theta(i,j) < 67.5)
                orientation = 2; 
            elseif (theta(i,j) >= 67.5) && (theta(i,j) <= 157.5)
                orientation = 3; 
            end
            switch orientation
                case 0 
                    if (mag(i,j) >= mag(i,j-1)) && (mag(i,j) >= mag(i,j+1))
                        edges(i,j) = mag(i,j);
                    end
                case 1 % 
                    if (mag(i,j) >= mag(i-1,j-1)) && (mag(i,j) >= mag(i+1,j+1))
                        edges(i,j) = mag(i,j);
                    end
                case 2 
                    if (mag(i,j) >= mag(i-1,j+1)) && (mag(i,j) >= mag(i+1,j-1))
                        edges(i,j) = mag(i,j);
                    end
                case 3 
                    if (mag(i,j) >= mag(i-1,j)) && (mag(i,j) >= mag(i+1,j))
                        edges(i,j) = mag(i,j);
                    end
            end
        end
    end
end