function b_image = freehand_blend(sourcename,targetname)
% use the algirthm shown in the link below to blend
%http://cs.brown.edu/courses/cs129/results/proj2/taox/

% read the image in the format of double precision 
source = im2double(imread(sourcename));
target = im2double(imread(targetname));

%show the image 
imshow(source, []);
%let user freehand select region they like
freehand = imfreehand();
fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(freehand,fcn);
wait(freehand);
% create binary mask for the freehand selection region
mask = freehand.createMask();

% get non-zero rows and non-zero cols of mask
[nonzero_rows, nonzero_cols] = find(mask);
%find the min and max row/col
min_row = min(nonzero_rows);
max_row = max(nonzero_rows);
min_col = min(nonzero_cols);
max_col = max(nonzero_cols);

%caculate the weight and height
w = max_col - min_col;
h = max_row - min_row;

%crop the mask/region from source by bounding box
source_cropped = source(min_row:max_row, min_col:max_col, :);
mask_cropped = mask(min_row:max_row, min_col:max_col);


% select the region in the target
figure; 
imshow(target);
rect = imrect(gca, [min_row min_col w h]);
setFixedAspectRatioMode(rect,1);
wait(rect);
des = rect.getPosition()

% get the size and offset of user-selected region in target image. 
xmin = des(1);
ymin = des(2);
rect_cols = des(3);
rect_rows = des(4);

%esize the cropped source and mask 
%and pad them to be the same dimensions as the target image
rf = imresize(source_cropped, [rect_rows, rect_cols]);
rm = imresize(mask_cropped, [rect_rows, rect_cols]);
source_padded = padarray(rf, [(size(target,1)-size(rf,1)), size(target,2)-size(rf,2)], 'post');
mask_padded = padarray(rm, [(size(target,1)-size(rm,1)), size(target,2)-size(rm,2)], 'post');

%translate the source and target rgarding to the offset
source = imtranslate(source_padded, [xmin,ymin]);
mask = imtranslate(mask_padded, [xmin,ymin]);

%find the new bouding box with its width and height
[nonzero_rows, nonzero_cols] = find(mask);
min_row = min(nonzero_rows);
max_row = max(nonzero_rows);
min_col = min(nonzero_cols);
max_col = max(nonzero_cols);
w = max_col - min_col;
h = max_row - min_row;

% corp the source, the mask, and the target based on the new bounding box
%all of them have the same width and height
source_region = source(min_row:max_row, min_col:max_col, :);
mask_region = mask(min_row:max_row, min_col:max_col);
target_region = target(min_row:max_row, min_col:max_col, :);

%% direct paste 
foreground = source .* repmat(mask, [1,1,3]);
background = target .* repmat(~mask, [1,1,3]);
direct_paste = background + foreground;
figure;
imshow(direct_paste);

%% blending
% use the algirthm shown in the link below
%http://cs.brown.edu/courses/cs129/results/proj2/taox/

%add a 1-pixel border around the source,target, and mask.
%will be removed after blending. 
source_region = padarray(source_region, [1,1], 'symmetric');
target_region = padarray(target_region, [1,1], 'symmetric');
mask_region = padarray(mask_region, [1,1]);

%get the size of target region
[rows, cols, ~] = size(target_region);

s = reshape(source_region, rows*cols, []);
t = reshape(target_region, rows*cols, []);

% initialize vector b
b = zeros(rows*cols, 3);
%row_vec has entries that represent row indexes of A
row_vec = zeros(rows*cols, 1);
%col_vec has entries that represent column indexes of A
col_vec = zeros(rows*cols, 1);
%value_vec stores values at specific positions inside A. 
value_vec = zeros(rows*cols, 1);

%track the current row sparse matrix A
eq_n = 1;

% iterate through rows of matrix A, 
%insert the appropriate values in the matrix A and b
for i = 1:rows*cols
    % If the current pixel location is in the mask
    if mask_region(i)
        b(i,:) = 4*s(i,:) - s(i-1,:) - s(i+1,:) - s(i+rows,:) - s(i-rows,:);
        
        % Insert a 4 into A at the i of the current central pixel.
        row_vec(eq_n) = i;
        col_vec(eq_n) = i;
        value_vec(eq_n) = 4;
        eq_n = eq_n + 1;
        
        %Insert a -1 for the pixel to four neighbours of the current pixel 
        row_vec(eq_n) = i;
        col_vec(eq_n) = i + 1;
        value_vec(eq_n) = -1;
        eq_n = eq_n + 1;

        row_vec(eq_n) = i;
        col_vec(eq_n) = i - 1;
        value_vec(eq_n) = -1;
        eq_n = eq_n + 1;

        row_vec(eq_n) = i;
        col_vec(eq_n) = i - rows;
        value_vec(eq_n) = -1;
        eq_n = eq_n + 1;
           
        row_vec(eq_n) = i;
        col_vec(eq_n) = i + rows;
        value_vec(eq_n) = -1;
        eq_n = eq_n + 1;
    else
        % If the current pixel location is not in the mask
        row_vec(eq_n) = i;
        col_vec(eq_n) = i;
        value_vec(eq_n) = 1;
        eq_n = eq_n + 1;
        
        b(i,:) = t(i,:);
    end
end

%create the sparse matrix 
A = sparse(row_vec, col_vec, value_vec, rows*cols, rows*cols);

toc
disp('Finished constructing the matrix A...')

% caculate the result for each channel
c_red = A \ b(:,1);
c_green = A \ b(:,2);
c_blue = A \ b(:,3);

% reshape to the original size
c_red = reshape(c_red, [rows, cols]);
c_green = reshape(c_green, [rows, cols]);
c_blue = reshape(c_blue, [rows, cols]);

% stack the channels back
result = zeros(rows, cols, 3);
result(:,:,1) = c_red;
result(:,:,2) = c_green;
result(:,:,3) = c_blue;

% remove the border:
result = result(2:rows-1, 2:cols-1, :);
%% show the final blending result
b_image = target;
b_image(min_row:max_row, min_col:max_col, :) = result(:,:,:);

figure; 
imshow(b_image);