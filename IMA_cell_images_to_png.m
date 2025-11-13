function IMA_cell_images_to_png(directory)
    % CONVERT_CELL_IMAGES_TO_PNG Convert Reference_Cell_*.mat files to PNG
    %   convert_cell_images_to_png(directory) loads all Reference_Cell_*.mat
    %   files from the specified directory, extracts RefStruct.Image,
    %   normalizes the images, and saves them as brightfield_Cell_*.png
    %
    %   Example:
    %       convert_cell_images_to_png('/path/to/data')
    %       convert_cell_images_to_png('.')  % current directory
    
    % Validate input
    if ~exist(directory, 'dir')
        error('Directory does not exist: %s', directory);
    end
    
    % Get all matching .mat files in the specified directory
    files = dir(fullfile(directory, 'Reference_Cell_*.mat'));
    
    if isempty(files)
        warning('No Reference_Cell_*.mat files found in: %s', directory);
        return;
    end
    
    % Loop through each file
    for i = 1:length(files)
        % Get the filename
        filename = files(i).name;
        full_path = fullfile(directory, filename);
        
        % Extract the cell number (xx part)
        tokens = regexp(filename, 'Reference_Cell_(.+)\.mat', 'tokens');
        cell_id = tokens{1}{1};
        
        % Load the .mat file
        data = load(full_path);
        img = data.RefStruct.Image;
        
        % Normalize the image
        img_normalized = (img - min(img(:))) / (max(img(:)) - min(img(:)));
        
        % Create output filename
        output_filename = fullfile(directory, sprintf('brightfield_Cell_%s.png', cell_id));
        
        % Save as PNG
        imwrite(img_normalized, output_filename);
        
        fprintf('Processed: %s -> %s\n', filename, sprintf('brightfield_Cell_%s.png', cell_id));
    end
    
    fprintf('Done! Processed %d images.\n', length(files));
end