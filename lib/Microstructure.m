classdef Microstructure
    %MICROSTRUCTURE Edited for Shear Flow Creep
    %   Detailed explanation goes here
    
    properties
       
        %Working Directory
        HomeDirectory;
        
        %Analysis Directory
        AnalysisDirectory;
        
        %AEH directory
        AEHDirectory;
        
        %Figure and axes handles
        FigureHandle;
        AxesHandle;
        
        %Analysis Tab handle arrays
        PhaseNumberHandle
        PhaseNameHandle;
        PhaseModeHandle;
        PhaseColorHandle;  
        PhaseNValueHandle;
        PhaseAdValueHandle;
        PhaseQValueHandle;
        
        %microsctructure color info
        Colors = {[1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [1 0 1];...
                  [0 1 1]; [1 165/255 0]; [190/255 190/255 190/255]};
        red = [1 0 0];
        green = [0 1 0];
        blue = [0 0 1];
        yellow = [1 1 0];
        magenta = [1 0 1];
        cyan = [0 1 1];
        orange = [1 165/255 0];
        gray = [190/255 190/255 190/255];
        
        %initial file info
        Filename;
        FileType;
        AutoSaveDirectory;
        ImageData;
        
        %image data
        DisplayData;
        DisplayDataIndices;
        NumberDataPoints;
        NumberPhases;
        
        %data point properties
        EBSDCorrectionMatrix;
        DataCoordinateList;
        DataEulerAngle;
        DataPhase;
        
        %phase properties
        PhaseName; 
        PhaseColorValue;
        PhaseAdValue;
        PhaseNValue;
        PhaseQValue;
        PhaseVolumeFraction;
        
        %analysis properties
        MeshDensity;
        CustomMeshDensityValue;
        InitialStrainRate;
        FinalStrainRate;
        NumberStrainRateIntervals;
        StrainRateIntervals;
        InitialTemperature;
        FinalTemperature;
        NumberTemperatureIntervals;
        DefaultStrainRateRange=0;
        InitialStrainRateTensor;
        FinalStrainRateTensor;
        DefaultEquiviscousValue;
        ShearTypeString;
        StrainRateVectors;
        LoadType=1;
        LastLoadType;
        PresetLoadNumber;
        PlotNewWindow;
       
        %FE information
        ThreeNodeCoordinateList;
        ThreeNodeElementIndexList;
        LastThreeNodeCoordinateList;
        LastThreeNodeElementIndexList;
        SixNodeCoordinateList;
        SixNodeElementIndexList;
        BoundaryNodeRelationsList;
        DataPointsInElementList;
        NumberElements;
        NumberNodes;
        ElementPhases;
        ElementGrains;
        LastElementGrains;
        NumberGrains;
        Grains;
        GrainsNormalized;
        GrainsReduced;
        GrainsMeshed;
        LastGrainsReduced;
        GrainPhases;  
        OriginalElementGrains;
        AllGrainEdges;
        AllGrainNodes;
        
        %Result Storage
        FieldString={'stress 11';'stress 22';'stress 33';'stress 32';'stress 31';'stress 12';...
            'max princ stress';'min princ stress';'max shear stress';'dev stress 11';'dev stress 22';'dev stress 33';...
            'max princ dev stress';'min princ dev stress';'2nd inv stress';...
            'strain rate 11';'strain rate 22';'strain rate 33';'strain rate 32';'strain rate 31';'strain rate 12';...
            'max princ strain rate';'min princ strain rate';'dev strain rate 11';'dev strain rate 22';'dev strain rate 33';...
            'max princ dev strain rate';'min princ dev strain rate';'2nd inv strain rate';...
            'viscosity';'power dissipation density'};
        OverallEffectiveAdValue;
        OverallEffectiveNValue;
        OverallEffectiveQValue;
        TemperatureAnalysisResults;
        TemperatureEffectiveAdValue;
        TemperatureEffectiveBValue;
        TemperatureEffectiveNValue;
        MacroShearStrainRate;
        TemperatureIntervals;
        MacroStressStrainRateData;
        MicroFieldData;
        MicroCoordinates;
        MicroStress;
        MicroStrain;
        MicroPowerDissipationDensity;
        MicroPrincStressDir;
        MicroPrincStrainDir;
        MicroViscosity;
        MacroStressInvariant;
        MacroStrainInvariant;
        MacroTemperature;
        MacroStress;
        MacroStrain;
        OverallR2;
        TemperatureR2;
        MacroscaleStressStrainRate
        
        %Radial Basis Function Data
        RBFGrainAMatrix;
        RBFGrainLMatrix;
        RBFGrainUMatrix;
        RBFGrainPMatrix;
        RBFGrainXData;
        RBFGrainYData;
        
        %Mesh Edit Variables
        EditFigureHandle
        EditAxesHandle;
        TriPlotHandle;
        NodePlotHandle;
        DeletedNodes;
        
        %New meshing items
        HoleGrains;
        GrainHoles;
        FineGrainPolylines;
        GrainPolylines;
        GrainJunctionPoints;
        
    end
    
    methods
        
    end
    
end

