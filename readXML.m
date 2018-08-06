clear all;
close all;
path = '/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/LGE/LGEcontour.cvi42wsx';
%cmr_file = xmlread(path);
con = CMR42ContourReader(path);
%% Read Dicom Image
[volume_image, slice_data, image_meta_data] = dicom23D();
%%
Im = zeros(size(volume_image,1), size(volume_image, 2));

% Epicardium
for c=1:length(con.contours(1).pts{2})
Im([(con.contours(1).pts{2}(c,1)*4)],[(con.contours(1).pts{2}(c,2)*4)])=1;
end

% Endocardium
for c=1:length(con.contours(1).pts{3})
Im([(con.contours(1).pts{3}(c,1)*4)],[(con.contours(1).pts{3}(c,2)*4)])=2;
end

figure();
imagesc(Im)

%% Without multiplication by 4 (floor)
Im_no4 = zeros(size(volume_image,1), size(volume_image, 2));

% Epicardium
for c=1:length(con.contours(1).pts{2})
Im_no4([floor(con.contours(1).pts{2}(c,1))],[floor(con.contours(1).pts{2}(c,2))])=1;
end

% Endocardium
for c=1:length(con.contours(1).pts{3})
Im_no4([floor(con.contours(1).pts{3}(c,1))],[floor(con.contours(1).pts{3}(c,2))])=2;
end

figure();
imagesc(Im_no4)
%% Look at original image
Img = zeros(1024, 1024);
Img(volume_image(:,:,1)>0)=1;

figure();
n = ceil(sqrt(size(volume_image, 3)));
for i = 1: size(volume_image, 3)
    
    subplot(n, n, i)
    imagesc(volume_image(:,:,i))
    colormap gray;
    title(strcat("Slice = ", num2str(i)))
end


%%
Im2 = zeros(1024, 1024);

% Epicardium
for c=1:length(con.contours(14).pts{2})
Im2([(con.contours(14).pts{2}(c,1)*4)],[(con.contours(14).pts{2}(c,2)*4)])=1;
end

% Endocardium
for c=1:length(con.contours(14).pts{3})
Im2([(con.contours(14).pts{3}(c,1)*4)],[(con.contours(14).pts{3}(c,2)*4)])=2;
end

figure();
imagesc(Im2)

%%
num_con = length(con.contours);
epi_flow = cell(1, num_con);
endo_flow = cell(1, num_con);
for i = 1:num_con
    ctype = con.contours(i).ctype;
    epi_flow{i} = zeros(1024);
    endo_flow{i} = zeros(1024);
    for j = 1:length(ctype)
        contour_type = ctype(j);
        if strcmp(contour_type{1}, 'saendocardialContour')
            for c = 1: length(con.contours(i).pts{j})
                endo_flow{i}([(con.contours(i).pts{j}(c,1)*4)], [(con.contours(i).pts{j}(c,2)*4)]) = 1;
            end
        elseif strcmp(contour_type{1}, 'saepicardialContour')
            for c = 1: length(con.contours(i).pts{j})
                epi_flow{i}([(con.contours(i).pts{j}(c,1)*4)], [(con.contours(i).pts{j}(c,2)*4)]) = 1;
            end
        end
    end    
end

%% Remove zero matrix
count = 1;
for i = 1 : num_con
    if any(epi_flow{i}(:))
        epi_flow_edit(:,:,count) = epi_flow{i}(:,:);
        endo_flow_edit(:,:,count) = endo_flow{i}(:,:);
        count = count + 1;
    end
end

%% Combine two contours
%contours = epi_flow_edit + endo_flow_edit;
%num = size(contours, 3);
%SE = strel('disk', 17);

%figure();
%for i = 1: num
%    
%    n = ceil(sqrt(num));
%    subplot(n, n, i)
%    closed_im = imclose(contours(:,:,i),SE);
%    imagesc(closed_im)
%end
%% Get center coordinate
contours = epi_flow_edit + endo_flow_edit;
num = size(contours, 3);
myocardium = zeros(1024, 1024, num);
center_array = zeros(num, 2);
for i = 1: num
    [x_epi, y_epi] = find(epi_flow_edit(:,:,i) ~= 0);
    x_epi_center = round(mean(x_epi));
    y_epi_center = round(mean(y_epi));

    epi_J = regiongrowing(epi_flow_edit(:,:,i),x_epi_center,y_epi_center);

    [x_endo, y_endo] = find(endo_flow_edit(:,:,i) ~= 0);
    x_endo_center = round(mean(x_endo));
    y_endo_center = round(mean(y_endo));

    endo_J = regiongrowing(endo_flow_edit(:,:,i),x_endo_center,y_endo_center);

    myocardium(:,:,i) = xor(epi_J, endo_J);
    center_array(i, :) = [(x_epi_center + x_endo_center)/2, (y_epi_center + y_endo_center)/2];
end

%% Show image
figure();
for i = 1:num
    subplot(n,n,i);
    imagesc(myocardium(:,:,i))
end



%%
figure();
z = 1:1:num;
scatter(center_array(:,1), center_array(:,2))
grid on;

%%
test_path = '/Users/jameszhang/Documents/Rohan/cvi42/circle_parser/james.cvi42wsx';
rect = CMR42ContourReader(test_path);

%% try rect contours
rect_coord1 = rect.contours(1).pts{1};
figure();
scatter(rect_coord1(:,1), rect_coord1(:,2))

rect_coord2 = rect.contours(2).pts{1};
rect_coord3 = rect.contours(2).pts{2};

figure();
scatter(rect_coord2(:,1), rect_coord2(:,2))
hold on;
scatter(rect_coord3(:,1), rect_coord3(:,2))

%% 2nd try
test_path1 = '/Users/jameszhang/Documents/Rohan/cvi42/circle_parser/1.cvi42wsx';
test_path2 = '/Users/jameszhang/Documents/Rohan/cvi42/circle_parser/2.cvi42wsx';
rect1 = CMR42ContourReader(test_path1);
rect2 = CMR42ContourReader(test_path2);

%%
coord1 = rect1.contours(1).pts{1};
coord2 = rect1.contours(2).pts{1};

figure();
scatter(coord1(:,1), coord1(:,2))
hold on;
scatter(coord2(:,1), coord2(:,2))

coord3 = rect2.contours(1).pts{1};
figure();
scatter(coord3(:,1), coord3(:,2))

%% 3rd try
test_path3 = '/Users/jameszhang/Documents/Rohan/cvi42/circle_parser/3.cvi42wsx';
rect3 = CMR42ContourReader(test_path3);
coord3 = rect3.contours(1).pts{1};
coord4 = rect3.contours(2).pts{1};

figure();
scatter(coord3(:,1), coord3(:,2))
hold on;
scatter(coord4(:,1), coord4(:,2))

%% 4th try
test_path4 = '/Users/jameszhang/Documents/Rohan/cvi42/circle_parser/4.cvi42wsx';
rect4 = CMR42ContourReader(test_path4);
