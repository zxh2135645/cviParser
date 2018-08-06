% second try
% From this try, we conclude that the XML is scaled in a way or shifted by
% some uncommon matrix size.
clear all;
close all;
path = '/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/LGE/LGEcontour.cvi42wsx';
%cmr_file = xmlread(path);
con = CMR42ContourReader(path);
%% Read Dicom Image
[volume_image, slice_data, image_meta_data] = dicom23D();
num_con = length(con.contours);
epi_flow = cell(1, num_con);
endo_flow = cell(1, num_con);
for i = 1:num_con
    ctype = con.contours(i).ctype;
    epi_flow{i} = zeros(size(volume_image, 1), size(volume_image,2));
    endo_flow{i} = zeros(size(volume_image, 1), size(volume_image,2));
    for j = 1:length(ctype)
        contour_type = ctype(j);
        if strcmp(contour_type{1}, 'saendocardialContour')
            for c = 1: length(con.contours(i).pts{j})
                endo_flow{i}([ceil(con.contours(i).pts{j}(c,1))], [ceil(con.contours(i).pts{j}(c,2))]) = 1;
            end
        elseif strcmp(contour_type{1}, 'saepicardialContour')
            for c = 1: length(con.contours(i).pts{j})
                epi_flow{i}([ceil(con.contours(i).pts{j}(c,1))], [ceil(con.contours(i).pts{j}(c,2))]) = 1;
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

%% Get center coordinate
contours = epi_flow_edit + endo_flow_edit;
num = size(contours, 3);
myocardium = zeros(size(volume_image,1), size(volume_image,2), num);
heart = zeros(size(volume_image,1), size(volume_image,2), num);
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
    center_array(i, :) = [(x_epi_center + x_endo_center)/2, (y_epi_center + y_endo_center)/2];
end

%% GetDICOMInfo
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
[volume_image, slice_data, image_meta_data] = ...
    dicom23D('/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/LGE/DICOM/HAM_CHOON_SIK_062Y/series0092-Body/', dicom_fields);
%info = dicominfo('/Users/jameszhang/Documents/Rohan/contour exporting/HAM_CHOON_SIK/LGE/DICOM/HAM_CHOON_SIK_062Y/series0092-Body/img0001-107.151.dcm');
%info.MediaStorageSOPInstanceUID
%% Find the index that contours matches with the image
for i = 1:length(con.contours)
    if strcmp(slice_data(1).MediaStorageSOPInstanceUID, con.contours(i).iuid)
        fprintf('%d \n',i)
    end
end

%% Visualize Myocardium
figure();
n = ceil(sqrt(num));
for i = 1:num
    subplot(n,n,i);
    imagesc(myocardium(:,:,i))
end
%% Visualize Heart
figure();
for i = 1:num
    subplot(n,n,i)
    imagesc(heart(:,:,i));
end
%% Apply Myocaridum to Image
figure();
mask = zeros(size(volume_image,1), size(volume_image,2));

for i = 1:num
    subplot(n,n,i);
    mask = myocardium(:,:,i) .* volume_image(:,:,1);
    imagesc(mask);
    colormap gray;
end

figure();
imagesc(volume_image(:,:,1))

%% Apply heart to Image

for i = 1:num
    figure();
    %subplot(n,n,i);
    mask_heart = heart(:,:,i) .* volume_image(:,:,1);
    imagesc(mask_heart);
    colormap gray;
end
%% To second slice
mask_heart = zeros(size(volume_image,1), size(volume_image,2));
figure();
for i = 1:num
    
    subplot(n,n,i);
    mask_heart = heart(:,:,i) .* volume_image(:,:,8);
    imagesc(mask_heart);
    colormap gray;
    alpha(0.5)
end
