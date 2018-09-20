%% Read excludeContour
clear all;
close all;

sequence_label = {'LGE', 'T1'};
label = char(sequence_label(2));

xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/LEE_KWAN_JOON/', label, '/*contour.cvi42wsx'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/LEE_KWAN_JOON/', label, '/*/*'));
cvi42wsx = xml_glob{1};
dicom = dicom_glob{1};
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
num_contours = length(slice_data);

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
excludeContour = cell(1, num_contours);
excludeCtr_struct = struct;

for i = 1:num_contours
    ctype = con.contours(contour_idx(i)).ctype;
    epi_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
    endo_flow{i} = zeros(size(volume_image, 2), size(volume_image,1));
    excludeContour{i} = zeros(size(volume_image, 2), size(volume_image,1));
    for j = 1:length(ctype)
        contour_type = ctype{j};
        if strcmp(contour_type, 'saendocardialContour')
            for c = 1: length(con.contours(contour_idx(i)).pts{j})
                endo_flow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
            end
        elseif strcmp(contour_type, 'saepicardialContour')
            for c = 1: length(con.contours(contour_idx(i)).pts{j})
                epi_flow{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
            end
        elseif contains(contour_type, 'excludeEnhancementAreaContour')
            for c = 1: length(con.contours(contour_idx(i)).pts{j})
                excludeContour{i}([ceil(con.contours(contour_idx(i)).pts{j}(c,1))], [ceil(con.contours(contour_idx(i)).pts{j}(c,2))]) = 1;
            end
            excludeCtr_struct.(contour_type) = excludeContour;
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

% For excludeContour
clear excludeContour_edit


fname = fieldnames(excludeCtr_struct);

for i = 1 : length(fname)
    count = 1;
    clear excludeCtr_mat
    ctr_index_array = [];
    for j = 1 : num_contours
        if any(excludeCtr_struct.(fname{i}){j}(:))
            excludeCtr_mat(:,:,count) =  excludeCtr_struct.(fname{i}){j};
            ctr_index_array = [ctr_index_array; j];
            count = count + 1;
        end
    end
    excludeContour_edit(i, 1:2) = {fname{i}, excludeCtr_mat};
    excludeContour_edit(i, 3) = {ctr_index_array};
end

%%
% Get center coordinate and region growing
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

shifted_heart = zeros(size(volume_image));
shifted_myo = zeros(size(volume_image));


for i = 1:num
    shifted_heart(:,:,i) = flipud(rot90(heart(:,:,i),1));
    shifted_myo(:,:,i) = flipud(rot90(myocardium(:,:,i),1));
    mask_heart = shifted_heart(:,:,i) .* volume_image(:,:,index_array(i));
    mask_myocardium = shifted_myo(:,:,i) .* volume_image(:,:,index_array(i));
end

%% Do I need info from index_array?
excludeContour = struct;

for i = 1:length(fname)
    figure();
    excludeContour.(fname{i}) = cell(1, 2);
    n = ceil(sqrt(size(excludeContour_edit{i, 2}, 3)));
    for j = 1:size(excludeContour_edit{i, 2}, 3)
        excludeContour_edit{i, 2}(:,:,j) = imfill(excludeContour_edit{i, 2}(:,:,j), 'holes');
        excludeContour_edit{i, 4}(:,:,j) = flipud(rot90(excludeContour_edit{i, 2}(:,:,j),1));
        subplot(n,n,j)
        imagesc(excludeContour_edit{i, 4}(:,:,j))
        axis equal
    end
    excludeContour.(fname{i}){1} = excludeContour_edit{i, 4};
    excludeContour.(fname{i}){2} = excludeContour_edit{i, 3};
end