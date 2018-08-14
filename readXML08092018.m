% Forth try for unmatched subjects - take HAN_DOL_BOON as an example
% UID doesn't match because the image is not the one we contoured
% 21 subjects are PSIR
clear all;
close all;

sequence_label = {'LGE', 'T1', 'PSIR', 'MAG'};
label = char(sequence_label(1));

xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/*/', label, '/*contour.cvi42wsx'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/*/', label, '/*/*'));

dicom = char(dicom_glob(1));
cvi42wsx = char(xml_glob(1));

strings = strsplit(dicom,'\');
CurrentFolder = pwd;
OutputPath = GetFullPath(cat(2, CurrentFolder, '/../masked/'));
dstFolder = cat(2, OutputPath, char(strings(end-4)), '/', label);

% Read DICOM and contours
con = CMR42ContourReader(cvi42wsx);
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
    
[volume_image, slice_data, image_meta_data] = dicom23D(dicom, dicom_fields);

%% Find Epi- and Endo-cardium from the contours
num_contours = length(con.contours);

epi_flow = cell(1, num_contours);
endo_flow = cell(1, num_contours);

for i = 1:num_contours
    ctype = con.contours(i).ctype;
    epi_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
    endo_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
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
index_array = [];
for i = 1 : num_contours
    if any(epi_flow{i}(:))
        epi_flow_edit(:,:,count) = epi_flow{i}(:,:);
        endo_flow_edit(:,:,count) = endo_flow{i}(:,:);
        index_array = [index_array; i];
        count = count + 1;
    end
end


%% Get center coordinate and region growing
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

%% Show image

for i = 1:num
    figure();
    shifted_heart(:,:,i) = flipud(rot90(heart(:,:,i),1));
    mask_heart = shifted_heart(:,:,i) .* volume_image(:,:,7);
    imagesc(mask_heart);
    axis equal;
end

%% 
sequence_label = {'LGE', 'T1', 'PSIR', 'MAG'};
label = char(sequence_label(3));

dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/EightSubjects_Guan/*/', label, '/*/*'));
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
        'SeriesDescription',...
        };

for i = 1: length(dicom_glob)
    dicom = char(dicom_glob(i));
    strings = strsplit(dicom,'\');
    name = char(strings(end-4));
    [volume_image, slice_data, image_meta_data] = dicom23D(dicom, dicom_fields);
    fprintf('Name: %s \n', name)
    fprintf('SeriesDescription: %s \n', slice_data(1).SeriesDescription)
end

%%
label = char(sequence_label(2));
%dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/EightSubjects_Guan/*/', label, '/*/*'));
%dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/', label, '/*/*'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/Wrong_Sequence/*/', label, '/*/*'));
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
        'SeriesDescription',...
        };

for i = 1: length(dicom_glob)
    dicom = char(dicom_glob(i));
    strings = strsplit(dicom,'\');
    name = char(strings(end-4));
    [volume_image, slice_data, image_meta_data] = dicom23D(dicom, dicom_fields);
    fprintf('Name: %s \n', name)
    fprintf('SeriesDescription: %s \n', slice_data(1).SeriesDescription)
end

%% Kim Tae Il as second example
clear all;
close all;


sequence_label = {'LGE', 'T1', 'PSIR', 'MAG'};
label = char(sequence_label(2));

xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/Wrong_Sequence/*/', label, '/*contour.cvi42wsx'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/Wrong_Sequence/*/', label, '/*/*'));
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
        'SeriesDescription',...
        };
    
dicom = char(dicom_glob(5));
strings = strsplit(dicom,'\');
name = char(strings(end-4));
[volume_image, slice_data, image_meta_data] = dicom23D(dicom, dicom_fields);
fprintf('Name: %s \n', name)
fprintf('SeriesDescription: %s \n', slice_data(1).SeriesDescription)

cvi42wsx = char(xml_glob(5));
con = CMR42ContourReader(cvi42wsx);

%%
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/KIM_TAE_IL/', label, '/*/*'));
dicom = char(dicom_glob(1));
strings = strsplit(dicom,'\');
name = char(strings(end-4));
[volume_image, slice_data, image_meta_data] = dicom23D(dicom, dicom_fields);
fprintf('Name: %s \n', name)
fprintf('SeriesDescription: %s \n', slice_data(1).SeriesDescription)
