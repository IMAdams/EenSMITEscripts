Params = [];
Params.MeanLevel = [];
Params.VarLevel = [];
Params.CCDVar = [];
Params.Gain = [];
Params.CCDOffset = [];
Params.CameraObj = [];

Params.MeanLevel = ParamsLg.Params.MeanLevel(129:384,129:384,:);
Params.VarLevel = ParamsLg.Params.VarLevel(129:384,129:384,:);
Params.CCDVar = ParamsLg.Params.CCDVar(129:384,129:384);
Params.Gain = ParamsLg.Params.Gain(129:384,129:384);
Params.CCDOffset = ParamsLg.Params.CCDOffset(129:384,129:384);
Params.CameraObj = ParamsLg.Params.CameraObj;

save("C:\Users\imadams\Documents\smite workspace\iX83_gainoffset\CH1_GainCal_256.mat", "Params");
