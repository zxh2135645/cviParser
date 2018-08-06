% Third try
% We found that we can match each contours to its original image, enabling
% us to move on.
% The coordinates are now rotated 90 degrees counter-clockwise 
% And then flip it upside down
% Succeed!

clear all;
close all;

% Read files
path = '/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/LGE/LGEcontour.cvi42wsx';
path = '/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/T1/T1contour.cvi42wsx';
con = CMR42ContourReader(path);
dicom_fields = {...
        'Filename',...
        'Height', ...
        'Width', ...
        'Rows',...
        'Columns', ...
        'PixelSpacing',...
        'SliceThickness',...
        'SliceLocation',...
        'SpacingBetweenSlices'...
        'ImagePositionPatient',...
        'ImageOrientationPatient',...
        'MediaStorageSOPInstanceUID',...
        };
%[volume_image, slice_data, image_meta_data] = ...
%    dicom23D('/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/LGE/DICOM/HAM_CHOON_SIK_062Y/series0092-Body/', dicom_fields);
[volume_image, slice_data, image_meta_data] = ...
    dicom23D('/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/T1/DICOM/HAM_CHOON_SIK_062Y/', dicom_fields);
%% Find Epi and Endo
num_contours = length(slice_data);
%num_contours = 6;

contour_idx = zeros(1, num_contours);

for i = 1:num_contours
    for j = 1:length(con.contours)
        if strcmp(slice_data(i).MediaStorageSOPInstanceUID, con.contours(j).iuid)
            contour_idx(i) = j;
        end
    end
end

epi_flow = cell(1, num_contours);
endo_flow = cell(1, num_contours);

for i = 1:num_contours
    ctype = con.contours(contour_idx(i)).ctype;
    epi_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
    endo_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
    for j = 1:length(ctype)
        contour_type = ctype(j);
        if strcmp(contour_type{1}, 'saendocardialContour')
            for c = 1: length(con.contours(contour_idx(i)).pts{j})
                endo_flow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
            end
        elseif strcmp(contour_type{1}, 'saepicardialContour')
            for c = 1: length(con.contours(contour_idx(i)).pts{j})
                epi_flow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
            end
        end
    end    
end

%% Remove zero matrix
clear endo_flow_edit epi_flow_edit;
count = 1;
for i = 1 : num_contours
    if any(epi_flow{i}(:))
        epi_flow_edit(:,:,count) = epi_flow{i}(:,:);
        endo_flow_edit(:,:,count) = endo_flow{i}(:,:);
        count = count + 1;
    end
end

%% Convert to matrix
%epi_flow_edit = zeros(size(volume_image, 2), size(volume_image, 1), num_contours);
%endo_flow_edit = zeros(size(volume_image, 2), size(volume_image, 1), num_contours);

%for i = 1:num_contours
%    epi_flow_edit(:,:,i) = epi_flow{i}(:,:);
%    endo_flow_edit(:,:,i) = endo_flow{i}(:,:);
%end

% Get center coordinate
contours = epi_flow_edit + endo_flow_edit;
num = size(contours, 3);
myocardium = zeros(size(volume_image,2), size(volume_image,1), num);
heart = zeros(size(volume_image,2), size(volume_image,1), num);
blood_pool = zeros(size(volume_image,2), size(volume_image,1), num);
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
    heart(:,:,i) = epi_J;
    blood_pool(:,:,i) = endo_J;
    center_array(i, :) = [(x_epi_center + x_endo_center)/2, (y_epi_center + y_endo_center)/2];
end

%{
%% Visualize Heart
n = ceil(sqrt(num));
%figure();
for i = 1:num
    %subplot(n,n,i)
    figure();
    imagesc(volume_image(:,:,i));
    colormap gray;
end
%%
for i = 1:num
    %subplot(n,n,i)
    figure();
    imagesc(blood_pool(:,:,i));
    colormap gray;
    axis equal
end
%}
%% Apply heart to Image
shifted_heart = zeros(size(volume_image));

for i = 1:num
    figure();
    %subplot(n,n,i);
    
    shifted_heart(:,:,i) = flipud(rot90(circshift(heart(:,:,i), [0, -0]),1));
    mask_heart = shifted_heart(:,:,i) .* volume_image(:,:,i);
    imagesc(mask_heart);
    %colormap gray;
    axis equal;
end

%% Save files to ../masked/
currentFolder = pwd;
dstFolder = GetFullPath(cat(2, currentFolder, '/../masked/'));

% need to import GetFullPath.m function

strings = strsplit(path,'/');
dstFolder = cat(2, dstFolder, char(strings(end-2)), '/', char(strings(end-1)));

if ~ exist(dstFolder, 'dir')
    mkdir(dstFolder);    
end

for i = 1:num
    mask_heart = shifted_heart(:,:,i) .* volume_image(:,:,i);
    dstPath = cat(2, dstFolder, '/masked_heart', num2str(i), '.mat');
    save(dstPath, 'mask_heart');
end