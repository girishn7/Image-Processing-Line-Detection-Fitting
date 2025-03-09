A = imread('C:\Users\giris\Downloads\hinges.pnm');
num_iterations = 1000;
inlier_threshold = 3;              
num_inliers_threshold = 50;        
[rows, cols] = find(A);     
best_line = [0, 0, 0];             
best_inliers = [];
for iter = 1:num_iterations
    rand_indices = randperm(length(rows), 2);
    point1 = [cols(rand_indices(1)), rows(rand_indices(1))];
    point2 = [cols(rand_indices(2)), rows(rand_indices(2))];
    slope = (point2(2) - point1(2)) / (point2(1) - point1(1));
    y_intercept = point1(2) - slope * point1(1);
    distances = abs(slope * cols - rows + y_intercept) / sqrt(slope^2 + 1);
    inlier_indices = find(distances < inlier_threshold);
    num_inliers = length(inlier_indices);
    if num_inliers > num_inliers_threshold && num_inliers > best_line(3)
        best_line = [slope, y_intercept, num_inliers];
        best_inliers = inlier_indices;
    end
end
x = linspace(1, size(A, 2), 100);
y = best_line(1) * x + best_line(2);
figure;
imshow(A);
title('RANSAC line detection')
hold on;
plot(cols(best_inliers), rows(best_inliers), 'g.');  
plot(x, y, 'r', 'LineWidth', 2);                    
hold off;
