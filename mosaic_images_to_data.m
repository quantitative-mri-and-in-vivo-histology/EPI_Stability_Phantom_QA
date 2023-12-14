function images = mosaic_images_to_data(dicom_list)
%MOSAIC_IMAGES_TO_DATA Summary of this function goes here
%   Detailed explanation goes here

for list_indx = 1:numel(dicom_list)
    image_info = dicominfo(dicom_list{list_indx});
    image_dcm = double(dicomread(dicom_list{list_indx}));
    
    % Mosaic Dimensionen ********Nur wenn phase H/F?********
    Mosaic_rows = image_info.Height / image_info.AcquisitionMatrix(4);
    Mosaic_columns = image_info.Width / image_info.AcquisitionMatrix(1);
    
    if ~isfield(image_info,'Private_0019_100a')
        slices = Mosaic_columns * Mosaic_rows;
    else
        slices = double(image_info.Private_0019_100a);
    end

    for slice_indx = 1:slices
        imgX = fix(double(slice_indx -1) / double(Mosaic_columns)) * image_info.AcquisitionMatrix(4) + 1;
        imgY = mod(slice_indx - 1,Mosaic_columns) * image_info.AcquisitionMatrix(1) +1;
        
        images(:,:,slice_indx,list_indx) = image_dcm(imgX:imgX+image_info.AcquisitionMatrix(4)-1,...
                                                     imgY:imgY+image_info.AcquisitionMatrix(1)-1);                                    
    end
    disp(['Loaded mosaic for repetition #' num2str(list_indx)]);
end

