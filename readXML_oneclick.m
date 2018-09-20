clear all;
close all;

sequence_label = {'LGE', 'T1'};
label = char(sequence_label(2));

xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/', '/contour.cvi42wsx'));
dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/contour_exporting_Guan/*/', label, '/*/*'));
%xml_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/CVIContour_HumanStudy/*/', label, '/*.cvi42wsx'));
%dicom_glob = glob(cat(2, 'C:/Users/ZhangX1/Desktop/CVIContour_HumanStudy/*/', label, '/DICOMS'));
CurrentFolder = pwd;
OutputPath = GetFullPath(cat(2, CurrentFolder, '/../masked/'));

%{ 
% Testing Eric's data 
dicom = char(dicom_glob(1));
cvi42wsx = char(xml_glob(1));
strings = strsplit(dicom,'\');
dstFolder = cat(2, OutputPath, char(strings(end-4)), '/', label);
con = CMR42ContourReader(cvi42wsx);
%}

dsts = {'Heart', 'Myocardium', 'excludeContour', 'MyoReference', 'MI'};
if length(xml_glob) ~= length(dicom_glob)
    errordlg('Number of DICOM does not match the number of contours, please check!')
else    
    for i = 1:length(xml_glob)
        dicom = char(dicom_glob(i));
        cvi42wsx = char(xml_glob(i));
        for dst = 1:length(dsts)
            strings = strsplit(dicom,'\');
            dstFolder = cat(2, OutputPath, char(strings(end-4)), '/', label, '/', dsts{dst});
            if ~ exist(dstFolder, 'dir')
                CMR42ContourMatrixGenerator(cvi42wsx, dicom, dstFolder, dsts{dst}, label)
            end
        end
    end
end

disp("Done!")
