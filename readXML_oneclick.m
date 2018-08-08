clear all;
close all;

sequence_label = {'LGE', 'T1'};
label = char(sequence_label(2));

xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/', label, '/*contour.cvi42wsx'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/', label, '/*/*'));
CurrentFolder = pwd;
OutputPath = GetFullPath(cat(2, CurrentFolder, '/../masked/'));

if length(xml_glob) ~= length(dicom_glob)
    errordlg('Number of DICOM does not match the number of contours, please check!')
else
    for i = 1:length(xml_glob)
        dicom = char(dicom_glob(i));
        cvi42wsx = char(xml_glob(i));
        strings = strsplit(dicom,'\');
        dstFolder = cat(2, OutputPath, char(strings(end-4)), '/', label);
        if ~ exist(dstFolder, 'dir')
            CMR42ContourMatrixGenerator(cvi42wsx, dicom, dstFolder)
        end
    end
end

