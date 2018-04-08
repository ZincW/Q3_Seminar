CourseworkC = struct;
CourseworkC.simg1 = im2double(imread('source1.jpg'));
CourseworkC_image = im2double(imread('target1.jpg'));
CourseworkC_cell = struct2cell(CourseworkC);
N=1;hold on
% if we need more than one image, just change the N
for n =1:N
figure(1)
CourseworkC_simg = CourseworkC_cell{n,1};
imshow(CourseworkC_simg,[]);
% if we want other selected region
% source = imfreehand();
source = imrect();
wait(source);
CourseworkC_scor = round(source.getPosition());
sxmin = CourseworkC_scor(1);
symin = CourseworkC_scor(2);
swidth = CourseworkC_scor(3);
sheight = CourseworkC_scor(4);
source = CourseworkC_simg(symin:(symin+sheight), sxmin:(sxmin+swidth), :);

%% crop the region of interest in destination image
% show the destination image
figure(2);
imshow(CourseworkC_image);
% des = imfreehand();
des = imrect();
wait(des);
dcor = round(des.getPosition());
drmin = dcor(1);
dcmin = dcor(2);
dwidth = dcor(3);
dheight = dcor(4);
resized_source = imresize(source,[dheight,dwidth]);
pad_source = padarray(resized_source, [(size(CourseworkC_image,1)-size(resized_source,1)),...
                                        size(CourseworkC_image,2)-size(resized_source,2)], 'post');
foreground = imtranslate(pad_source, [drmin,dcmin]);

[Sr,Sc] = size(foreground(:,:,1));
[r,c]= find(foreground(:,:,1));
rmin = min(r);
rmax = max(r);
cmin = min(c);
cmax = max(c);

% now we construct the background
background = CourseworkC_image;
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
Uimd = reshape(CourseworkC_image,Sr*Sc,3);
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
a = 1e10;
terml = G'*G + a*S'*S;
termr = G'*gtilda_update + a*S'*S*Uimd;
termr = sparse(termr);
% get final result from linear equation
[L,R,P] = lu(terml);
z = [];
for i = 1:3
    z(:,i) = L\(P*termr(:,i));  
end
U = [];
for i = 1:3
    U(:,i) = R\z(:,i);
end
CourseworkC_image = reshape(U,Sr,Sc,3);
print('finished!')
end
figure(4)
imshow(uint8(CourseworkC_image*255));
