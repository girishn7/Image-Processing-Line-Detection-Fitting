A = imread('C:\Users\giris\Downloads\hinges.pnm'); 
sigma = 2; 
low_thresh = 0.0001; 
high_thresh = 0.2; 
subplot(2,2,1), imshow(A), title('Original Image');
edges = canny_edge_detection(A, sigma, low_thresh, high_thresh); 

theta_resolution = 1; 
theta = -90:theta_resolution:89;
max_rho = round(sqrt(size(edges, 1)^2 + size(edges, 2)^2)); 
rho_resolution = 1; 
num_rhos = round(2*max_rho/rho_resolution) + 1;
H = zeros(num_rhos, length(theta));

[rows, cols] = find(edges);
num_edges = length(rows);

for i = 1:num_edges
    x = cols(i);
    y = rows(i);
    for j = 1:length(theta)
        rho = x*cosd(theta(j)) + y*sind(theta(j));
        rho_index = round(rho/rho_resolution) + max_rho/rho_resolution + 1;
        H(rho_index, j) = H(rho_index, j) + 1;
    end
end


subplot(2,2,3);
imagesc(theta, (-max_rho:max_rho)*rho_resolution, H);
colormap(hot);
colorbar;
xlabel('Theta (degrees)');
ylabel('Rho');
title('Hough Accumulator Array');

threshold = 0.3*max(H(:)); 
peaks = H > threshold;
[rho_indices, theta_indices] = find(peaks);

line_segments = cell(length(rho_indices), 1);
for i = 1:length(rho_indices)
    rho = (rho_indices(i) - max_rho/rho_resolution - 1)*rho_resolution;
    th = theta(theta_indices(i));
    x = 0:size(edges, 2)-1;
    y = round((rho - x*cosd(th))/sind(th));
    line_segments{i} = [x', y'];
end

line_parameters = zeros(length(line_segments), 2); 
for i = 1:length(line_segments)
    coords = line_segments{i};
    mean_coords = mean(coords, 1);
    centered_coords = coords - mean_coords;
    cov_matrix = centered_coords'*centered_coords;
    [eig_vectors, eig_values] = eig(cov_matrix);
    [~, ind] = min(diag(eig_values));
    line_dir = eig_vectors(:, ind);
    line_dir = line_dir/norm(line_dir);
    line_parameters(i, :) = [mean_coords*line_dir, atan2d(line_dir(2), line_dir(1))];
end

subplot(2,2,4);
imshow(A); title('Line Fitting');
hold on;
for i = 1:length(line_segments)
    segment = line_segments{i};
    plot(segment(:, 1), segment(:, 2), 'LineWidth', 2);
end


 
function edges = canny_edge_detection(A, sigma, low_thresh, high_thresh)


filter_size = 2*ceil(3*sigma)+1; 
filter = fspecial('gaussian', filter_size, sigma); 
A_smoothed = conv2(A, filter, 'same'); 

dx_filter = [-1 0 1; -sqrt(2) 0 sqrt(2); -1 0 1]/(2*sqrt(2)); 
dy_filter = dx_filter'; 
dx = conv2(A_smoothed, dx_filter, 'same'); 
dy = conv2(A_smoothed, dy_filter, 'same'); 
mag = sqrt(dx.^2 + dy.^2); 
theta = atan2d(dy, dx);

[m, n] = size(A);
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
            case 1 
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


edges_high = edges > high_thresh; 
edges_low = edges > low_thresh;
edges_low(1,:) = 0; edges_low(m,:) = 0; edges_low(:,1) = 0; edges_low(:,n) = 0; 
edges_strong = edges_high;
for i = 2:m-1
    for j = 2:n-1
        if edges_high(i,j) && ~edges_strong(i,j)
            
            neighbors = edges_low(i-1:i+1,j-1:j+1);
            edges_strong(i,j) = any(neighbors(:));
        end
    end
end
edges = edges_strong;
subplot(2,2,2), imshow(edges), title('Edge Map input for Hough');
end
