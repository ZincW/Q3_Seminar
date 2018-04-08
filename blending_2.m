%% crop the region of interest of source image
% make a struct for all the source images
file = struct;
file.simg1 = im2double(imread('source.jpeg'));
dimg = im2double(imread('target.jpg'));
cell = struct2cell(file);
% N is the number of input source images
N= 1;
hold on
%Loop for all the source images
for n =1:N
%crop the region of interest in source image
figure(1)
% imread the nth source image
simg = cell{n,1};
imshow(simg,[]);
% source = imfreehand();

[t_rows, t_cols, channels] = size(dimg);

% Open an interactive session to allow the user to freehand select a region
% of interest in the source image, and create a mask for that region:
mask = create_freehand_mask(simg);
%wait(mask);
% We find the bounding box containing the nonzero values in the mask:
[nonzero_rows, nonzero_cols] = find(mask);
min_row = min(nonzero_rows);
max_row = max(nonzero_rows);
min_col = min(nonzero_cols);
max_col = max(nonzero_cols);

width = max_col - min_col;
height = max_row - min_row;

% We crop the source image and the mask itself to the size of the bounding
% box, since the only pixels we care about from the source are those within
% the mask:
source_cropped = simg(min_row:max_row, min_col:max_col, :);
mask_cropped = mask(min_row:max_row, min_col:max_col);

figure; imshow(source_cropped);
figure; imshow(mask_cropped);

%% crop the region of interest in destination image
% show the destination image
figure(2);
imshow(dimg);

target = dimg

des = imrect(gca, [min_row min_col width height]);
setFixedAspectRatioMode(des,1);
wait(des);

% get the coordinates of detination region
dcor = round(des.getPosition());
drmin = dcor(1);
dcmin = dcor(2);
dwidth = dcor(3);
dheight = dcor(4);

xmin = drmin;
ymin = dcmin;
rect_cols = dwidth;
rect_rows = dheight;

% Resize the source and mask if desired:
rf = imresize(source_cropped, [rect_rows, rect_cols]);
rm = imresize(mask_cropped, [rect_rows, rect_cols]);

% Pad the source and the mask to be the same dimensions as the target image
source_padded = padarray(rf, [(size(target,1)-size(rf,1)), size(target,2)-size(rf,2)], 'post');
mask_padded = padarray(rm, [(size(target,1)-size(rm,1)), size(target,2)-size(rm,2)], 'post');


% Translate the source and target by the user-given (x,y) offset
source = imtranslate(source_padded, [xmin,ymin]);
mask = imtranslate(mask_padded, [xmin,ymin]);

% Now we find the new bounding box for the translated and resized region:
[nonzero_rows, nonzero_cols] = find(mask);
min_row = min(nonzero_rows);
max_row = max(nonzero_rows);
min_col = min(nonzero_cols);
max_col = max(nonzero_cols);

width = max_col - min_col;
height = max_row - min_row;

% Now we crop all the images: the source, the mask, and the target, to
% focus on the region that has to be blended:
source_region = source(min_row:max_row, min_col:max_col, :);
mask_region = mask(min_row:max_row, min_col:max_col);
target_region = target(min_row:max_row, min_col:max_col, :);

%%
foreground = source .* repmat(mask, [1,1,3]);
background = target .* repmat(~mask, [1,1,3]);
direct_copy = background + foreground;
figure;
imshow(direct_copy);
%%

% Translate the source to the corresponding location in destination image

[Sr,Sc] = size(foreground(:,:,1));
[r,c]= find(foreground(:,:,1));
rmin = min(r);
rmax = max(r);
cmin = min(c);
cmax = max(c);

% % make the background
% background = dimg;
% for i = rmin:rmax
%     for j = cmin:cmax
%         for rgb = 1:3
%         background(i,j,rgb) = 0;
%         end
%     end
% end
% 
% % simple copy of the two images
% combine = foreground + background;
% figure(3)
% imshow(combine);
%% Blend the images

% Calculate gradient matrix
Uim1 = reshape(foreground, Sr*Sc, 3);
Uimd = reshape(dimg,Sr*Sc,3);
G = gradient(Sr,Sc);
gtilda = G*Uim1;
% set boundary gradient to zero
gtilda_update = coordinate(foreground,gtilda);
% Calculate S matrix
S = ones(Sr,Sc);
for i = rmin:rmax
    for j = cmin:cmax
       S(i,j) = 0; 
    end
end
S = build_S(S);
% set coefficient a
a = 1e10;
terml = G'*G + a*S'*S;
termr = G'*gtilda_update + a*S'*S*Uimd;
termr = sparse(termr);
% solve linear equation
[L,R,P] = lu(terml);
z = [];
for i = 1:3
    z(:,i) = L\(P*termr(:,i));  
end
U = [];
for i = 1:3
    U(:,i) = R\z(:,i);
end

dimg = reshape(U,Sr,Sc,3);
print('finished!')
end
figure(4)
imshow(uint8(dimg*255));