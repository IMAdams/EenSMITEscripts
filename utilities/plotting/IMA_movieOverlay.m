ImageStack = zeros([size(GaussIm1, 1), size(GaussIm1, 2), 2]);
      ImageStack(:, :, 1) = GaussIm1;
      ImageStack(:, :, 2) = GaussIm2;
      [OverlayImage, ColorOrderTag] = ...
         smi_vis.GenerateImages.overlayNImages(ImageStack);
      imwrite(OverlayImage, hot(256), [SaveFile, ‘.png’]);