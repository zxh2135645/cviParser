function CMR42ContourMatrixGenerator(cvi42wsx, dicom, dstFolder)

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

% Find Epi- and Endo-cardium from the contours
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

% Remove zero matrix
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


% Apply heart to the image and save it to OutputPath

if ~ exist(dstFolder, 'dir')
    mkdir(dstFolder);    
end

for i = 1:num
    shifted_heart(:,:,i) = flipud(rot90(heart(:,:,i),1));
    mask_heart = shifted_heart(:,:,i) .* volume_image(:,:,index_array(i));
    dstPath = cat(2, dstFolder, '/masked_heart', num2str(index_array(i)), '.mat');
    save(dstPath, 'mask_heart');
end

end