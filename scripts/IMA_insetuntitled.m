%% make an inset

BaseName = "FC_on_DC_on_LoT";
PlotSaveDir1 = "C:\Users\imadams\Documents\GroundZero\ian\SR_lab\Results\insets";
[GaussIm] = smi_vis.GenerateImages.gaussianImage(SMD, 20);
GaussImInset = GaussIm(1200:2200,1800:3400,:);

GIMF = figure();
GIMFA = axes(GIMF);
imshow(GaussImInset, [], 'Parent', GIMFA)
FileName = strcat(BaseName, '_GaussINSET.png');
imwrite(GaussImInset, fullfile(PlotSaveDir1, FileName))