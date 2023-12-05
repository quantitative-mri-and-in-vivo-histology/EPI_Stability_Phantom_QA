function images = single_slice_to_data(dicom_list)

%-FJv23(12.04): "Enhanced DICOM" modification added (all volumes added at
%once).
slice_pos_list = [];
slice_pos_counter = [];
images = [];
reordering_needed = 0;

for list_indx = 1:numel(dicom_list)
    image_info = dicominfo(dicom_list{list_indx});
    image_dcm = dicomread(dicom_list{list_indx});
    
    if numel(size(squeeze(image_dcm))) == 3
        %Here it is assumed that the slices are all loaded AND in the
        %correct order.
        images(:,:,:,list_indx) = double(squeeze(image_dcm));
        disp(['Loaded volume for repetition #' num2str(list_indx)]);
        
    else %Here it is assumed that each dcm file is a slice/repetition
        reordering_needed = 1;
        slice_pos = image_info.SliceLocation;

        if isempty(slice_pos_list == slice_pos) || isempty(find(slice_pos_list == slice_pos))
            slice_pos_list = [slice_pos_list, slice_pos];
            slice_pos_counter = [slice_pos_counter, 1];
            images(:,:,numel(slice_pos_list),1) = double(image_dcm);
            counter_disp = 1;
        else
            slice_pos_indx = find(slice_pos_list == slice_pos);
            slice_pos_counter(1,slice_pos_indx) = slice_pos_counter(1,slice_pos_indx) + 1;

            images(:,:,slice_pos_indx,slice_pos_counter(1,slice_pos_indx)) = double(image_dcm);
            counter_disp = slice_pos_counter(1,slice_pos_indx);
        end
        disp(['Loaded image in slice position: ' num2str(slice_pos) ', repetition #' num2str(counter_disp)]);
    end
end

if reordering_needed
    [~,indx_slice] = sort(slice_pos_list);
    images = images(:,:,indx_slice,:);
end

end