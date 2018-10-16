clear all;
close all;

sequence_label = {'LGE', 'T1'};
label = char(sequence_label(2));

xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/James_ground_truth/', '*.cvi42wsx'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/', label, '/*/*'));
CurrentFolder = pwd;
OutputPath = GetFullPath(cat(2, CurrentFolder, '/../CNNData/'));

dsts = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI', 'noReflowAreaContour'};
Names = cell(length(xml_glob), 1);
for i = 1:length(xml_glob)
    strings = strsplit(xml_glob{i},'\');
    name_ext = strsplit(strings{end}, '.');
    name = name_ext{1};
    Names{i} = name;
end

FullNames = cell(length(dicom_glob), 1);
for i = 1:length(dicom_glob)
    strings = strsplit(dicom_glob{i},'\');
    name = strings{end-4};
    FullNames{i} = name;
end

Name_labels = zeros(length(FullNames), 1);
for i = 1:length(Names)
    idx = find(strcmp(Names{i}, FullNames));
    Name_labels(idx) = i;
end
%%
for i = 1:length(Name_labels)
   if Name_labels(i) ~= 0
       dicom = dicom_glob{i};
       cvi42wsx = xml_glob{Name_labels(i)};
       
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
       for dst = 1:length(dsts)
           dstFolder = cat(2, OutputPath, FullNames{i}, '\', label, '\', dsts{dst});
           CMR42ContourMatrixGenerator2(con, volume_image, slice_data, dstFolder, dsts(dst), label)
       end
   end
end


disp("Done!")