function writeTiff16(data,fname)
% This is used to save the tif images as floating point values. The builtin
% functions do not work properly and we must make use of the Tiff writing
% function.
dims=size(data);
if numel(dims)<3
    maxZ=1;
    flag=1;
else
    maxZ=dims(3);
    flag=0;
end
for zz=1:maxZ
    switch zz
        case 1
            opt='w';
        otherwise
            opt='a';
    end
    
    if flag
        A=double(data);
    else
        % This is used to save the tif images as floating point values.
        A=double(data(:,:,zz-1));
    end
    t = Tiff([fname,'.tif'], opt);
    tagstruct.ImageLength = size(A, 1);
    tagstruct.ImageWidth = size(A, 2);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.ResolutionUnit=Tiff.ResolutionUnit.Inch;
    tagstruct.BitsPerSample = 64;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tagstruct);
    t.write(A);
    t.close();
end

end