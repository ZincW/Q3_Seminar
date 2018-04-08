%% crop the region of interest of source image
% make a struct for all the source images
file = struct;
file.simg1 = im2double(imread('source.jpg'));
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

source = imrect();
wait(source);
% get the coordinates of the selected region
scor = round(source.getPosition());
sxmin = scor(1);
symin = scor(2);
swidth = scor(3);
sheight = scor(4);
source = simg(symin:(symin+sheight), sxmin:(sxmin+swidth), :);


%% crop the region of interest in destination image
% show the destination image
figure(2);
imshow(dimg);
% des = imfreehand();
des = imrect();
wait(des);
% get the coordinates of detination region
dcor = round(des.getPosition());
drmin = dcor(1);
dcmin = dcor(2);
dwidth = dcor(3);
dheight = dcor(4);

% resize the source region to the same size of the destination region
resized_source = imresize(source,[dheight,dwidth]);

% make the source and mask have the same dimension as the destination image
pad_source = padarray(resized_source, [(size(dimg,1)-size(resized_source,1)),...
                                        size(dimg,2)-size(resized_source,2)], 'post');
% Translate the source to the corresponding location in destination image
foreground = imtranslate(pad_source, [drmin,dcmin]);

[Sr,Sc] = size(foreground(:,:,1));
[r,c]= find(foreground(:,:,1));
rmin = min(r);
rmax = max(r);
cmin = min(c);
cmax = max(c);

% make the background
background = dimg;
for i = rmin:rmax
    for j = cmin:cmax
        for rgb = 1:3
        background(i,j,rgb) = 0;
        end
    end
end

% simple copy of the two images
combine = foreground + background;
figure(3)
imshow(combine);
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
