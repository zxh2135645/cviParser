function pat = CMR42ExtractPatientData(cmr_file)
% This function reads CMR42 XML file and extracts only the PatientData key
%
%   pat = CMR42ExtractPatientData(cmr_file);
%
% Author: Avan Suinesiaputra - University of Auckland 2016

try
    cmr42XML = xmlread(cmr_file);
catch
    error('Cannot read %s file as an XML file.',cmr_file);
end

pat = struct;

allItems = cmr42XML.getElementsByTagName('Hash:item');
for k=0:allItems.getLength-1
    
    if( strcmp(allItems.item(k).getAttribute('Hash:key'),'PatientData') )
        
        patientData = allItems.item(k);
        for i=0:patientData.getLength-1
            
            % check hash:item
            if( ~strcmp(patientData.item(i).getNodeName, 'Hash:item') )
                continue;
            end
            
            % the strings are in java.lang.string
            % convert to matlab's string with char() function
            
            attr = patientData.item(i).getAttribute('Hash:key');
            type = patientData.item(i).getAttribute('Variant:type');
            if( strcmp(type,'double') )
                pat.(char(attr)) = str2double(patientData.item(i).getTextContent);
            elseif( strcmp(type,'QString') )
                strdata = char(patientData.item(i).getTextContent);
                if( ~isempty(strtrim(strdata)) )
                    pat.(char(attr)) = strdata;
                end
            end
            
        end
        
        break;    
    end
    
end