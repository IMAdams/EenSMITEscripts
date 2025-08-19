%% Generate ROI overlays with clusters

close all

options = {'MAPN', 'Gaussian', 'Boundary', 'Cluster', 'OneImage'};
start_DataDir = 'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays';
PixelSize = 97.4;
SaveDir = 'C:\Users\imadams\Documents\smite workspace\ROI cluster overlays\wtRestOverlay';

CI = smi_cluster.ClusterInterface;
CI.plotROIDriver(PixelSize, options, start_DataDir, SaveDir);
