clc;    clear;  close all;
%% Model selection
name_model = 'Houston_AL_canesm_20thcal';       % CanESM, present
% name_model = 'Houston_AL_canesm_ssp585cal';     % CanESM, future
% name_model = 'Houston_AL_gfdl6_20thcal';        % GFDL, present
% name_model = 'Houston_AL_gfdl6_ssp585cal';      % GFDL, future
% name_model = 'Houston_AL_hadgem6_20thcal';      % HadGEM, present
% name_model = 'Houston_AL_hadgem6_ssp585cal';    % HadGEM, future

path_model = [pwd,'\',name_model];
load([path_model,'\',name_model,'.mat']);       % load dataset
load([path_model,'\vnet.mat']);                 % load velocity data

%% Model values
Manning_n = 0.015;          % Default Manning's n
p_n = 101500;               % background atmospheric pressure
MHHW = 0.18;                % Mean higher-high water level
GSLR = 0;                   % No global SLR (present)
% GSLR = 0.48;              % Low global SRL
% GSLR = 0.57;              % Middle global SRL
% GSLR = 0.75;              % High global SRL
RSLR = 0;                   % No regional SLR (present)
% RSLS = 0.51;              % Regional SLR
h0 = GSLR + RSLS + MHHW;    % Initial water level condition

% Air-sea drag coefficient (Zweers et al. 2010)
C1 = 0.0012;
C2 = 0.003;
C3 = 0.0015;
W1 = 0;
W2 = 28;
W3 = 60;

%% TC filtering
lim1 = 96*0.44704;
lim2 = 111*0.44704;
lim3 = 130*0.44704;
lim4 = 157*0.44704;
lon_limit = -93.87;
lat_limit = 29.55;
cat_limit = 27.5;
t_interval = 3600;
vnet_max = zeros(4500,1);
for j=1:4500
    k = find(latstore(j,:)>lat_limit,1);
    if isempty(k)==0
        if longstore(j,k)-360<lon_limit
        	k_lim = find(latstore(j,:)>cat_limit);
            vnet_max(j) = max(vnet(j,k_lim))*0.51444;
        end
    end
end
track_list = transpose(find(vnet_max>=lim2));

%% Model setup
for track_num = track_list(1:end)
    idx_start = find(latstore(track_num,:)>=22,1);          % start time point
    idx_final = find(monthstore(track_num,:)==0,1);         % end time point
    year = double(yearstore(1,track_num));
    month = monthstore(track_num,idx_start:idx_final-1);
    day = daystore(track_num,idx_start:idx_final-1);
    hour = hourstore(track_num,idx_start:idx_final-1);
    time = datenum(year,month,day,hour,0,0)*86400;
    
    path_year = [path_model,'\year_',num2str(year)];
    path_track = [path_year,'\track_',num2str(track_num)];
    mkdir(path_track);
    
    %% Write hydrodynamic model files
    path_flow = [path_track,'\dflowfm'];
    mkdir(path_flow);    
    copyfile([pwd,'\dflowfm\*'],path_flow);
    
    % Copy spiderweb file (.spw)
    name_spider = ['Wind_grid_',num2str(track_num),'.spw'];
    copyfile([path_model,'\',name_spider],path_flow);
    
    % Write Extension file (.ext)
    name_ext = [path_flow,'\FlowFM.ext'];
    h_ext = fopen(name_ext,'w');
    fprintf(h_ext,'\n');
    fprintf(h_ext,'QUANTITY=frictioncoefficient \n');
    fprintf(h_ext,'FILENAME=frictioncoefficient.xyz \n');
    fprintf(h_ext,'FILETYPE=7 \n');
    fprintf(h_ext,'METHOD=6 \n');
    fprintf(h_ext,'OPERAND=O \n');
    fprintf(h_ext,'AVERAGINGTYPE=2 \n');
    fprintf(h_ext,'RELATIVESEARCHCELLSIZE=1 \n\n');
    fprintf(h_ext,'QUANTITY=airpressure_windx_windy \n');
    fprintf(h_ext,['FILENAME=',name_spider,'\n']);
    fprintf(h_ext,'FILETYPE=5 \n');
    fprintf(h_ext,'METHOD=1 \n');
    fprintf(h_ext,'OPERAND=+ \n');
    fclose(h_ext);
    
    % Write input file (.mdu)
    name_mdu = [path_flow,'\FlowFM.mdu'];
    h_mdu = fopen(name_mdu,'w');
    fprintf(h_mdu,'[General]\n');
    fprintf(h_mdu,'Program = D-Flow FM\n');
    fprintf(h_mdu,'Version = 1.2.110.67911M\n');
    fprintf(h_mdu,'FileVersion = 1.02\n');
    fprintf(h_mdu,'GuiVersion = 4.3.0.0\n');
    fprintf(h_mdu,'AutoStart = 0\n');
    fprintf(h_mdu,'PathsRelativeToParent = 1\n\n');

    fprintf(h_mdu,'[geometry]\n');
    fprintf(h_mdu,'NetFile = FlowFM_net.nc\n');
    fprintf(h_mdu,'BathymetryFile = \n');
    fprintf(h_mdu,'DryPointsFile = \n');
    fprintf(h_mdu,'GridEnclosureFile = \n');
    fprintf(h_mdu,'WaterLevIniFile = \n');
    fprintf(h_mdu,'LandBoundaryFile = \n');
    fprintf(h_mdu,'ThinDamFile = \n');
    fprintf(h_mdu,'FixedWeirFile = FlowFM_fxw.pliz \n');
    fprintf(h_mdu,'PillarFile = \n');
    fprintf(h_mdu,'StructureFile = \n');
    fprintf(h_mdu,'VertplizFile = \n');
    fprintf(h_mdu,'ProflocFile = \n');
    fprintf(h_mdu,'ProfdefFile = \n');
    fprintf(h_mdu,'ProfdefxyzFile = \n');
    fprintf(h_mdu,'Uniformwidth1D = 2 \n');
    fprintf(h_mdu,'ManholeFile = \n');
    fprintf(h_mdu,'WaterLevIni = %1.4f \n',h0);
    fprintf(h_mdu,'Bedlevuni = -5 \n');
    fprintf(h_mdu,'Bedslope = 0 \n');
    fprintf(h_mdu,'BedlevType = 3 \n');
    fprintf(h_mdu,'Blmeanbelow = -999 \n');
    fprintf(h_mdu,'Blminabove = -999 \n');
    fprintf(h_mdu,'PartitionFile = \n');
    fprintf(h_mdu,'AngLat = 0 \n');
    fprintf(h_mdu,'AngLon = 0 \n');
    fprintf(h_mdu,'Conveyance2D = -1 \n');
    fprintf(h_mdu,'Nonlin2D = 0 \n');
    fprintf(h_mdu,'Sillheightmin = 0 \n');
    fprintf(h_mdu,'Makeorthocenters = 0 \n');
    fprintf(h_mdu,'Dcenterinside = 1 \n');
    fprintf(h_mdu,'Bamin = 1E-06 \n');
    fprintf(h_mdu,'OpenBoundaryTolerance = 3 \n');
    fprintf(h_mdu,'RenumberFlowNodes = 1 \n');
    fprintf(h_mdu,'Kmx = 0 \n');
    fprintf(h_mdu,'Numtopsig = 0 \n');
    fprintf(h_mdu,'SigmaGrowthFactor = 1.2 \n');
    fprintf(h_mdu,'UseCaching = 1 \n');
    fprintf(h_mdu,'DzTop = 1 \n');
    fprintf(h_mdu,'FloorLevTopLay = -1.0 \n');
    fprintf(h_mdu,'DzTopUniAboveZ = -5.0 \n');
    fprintf(h_mdu,'NumTopSigUniform = 0 \n\n');

    fprintf(h_mdu,'[numerics]\n');
    fprintf(h_mdu,'CFLMax = 0.7 \n');
    fprintf(h_mdu,'AdvecType = 33 \n');
    fprintf(h_mdu,'TimeStepType = 2 \n');
    fprintf(h_mdu,'Limtyphu = 0 \n');
    fprintf(h_mdu,'Limtypmom = 4 \n');
    fprintf(h_mdu,'Limtypsa = 4 \n');
    fprintf(h_mdu,'Icgsolver = 2 \n');
    fprintf(h_mdu,'Maxdegree = 6 \n');
    fprintf(h_mdu,'FixedWeirScheme = 9 \n');
    fprintf(h_mdu,'FixedWeirContraction = 1 \n');
    fprintf(h_mdu,'Izbndpos = 0 \n');
    fprintf(h_mdu,'Tlfsmo = 3600 \n');
    fprintf(h_mdu,'Slopedrop2D = 0 \n');
    fprintf(h_mdu,'Chkadvd = 0.1 \n');
    fprintf(h_mdu,'Teta0 = 0.55 \n');
    fprintf(h_mdu,'Qhrelax = 0.01 \n');
    fprintf(h_mdu,'cstbnd = 0 \n');
    fprintf(h_mdu,'Maxitverticalforestersal = 0 \n');
    fprintf(h_mdu,'Maxitverticalforestertem = 0 \n');
    fprintf(h_mdu,'Turbulencemodel = 3 \n');
    fprintf(h_mdu,'Turbulenceadvection = 3 \n');
    fprintf(h_mdu,'AntiCreep = 0 \n');
    fprintf(h_mdu,'Maxwaterleveldiff = 0 \n');
    fprintf(h_mdu,'Maxvelocitydiff = 0 \n');
    fprintf(h_mdu,'Epshu = 0.0001 \n');
    fprintf(h_mdu,'FixedWeirRelaxationcoef = 0.6 \n\n');

    fprintf(h_mdu,'[physics]\n');
    fprintf(h_mdu,'UnifFrictType = 1 \n');
    fprintf(h_mdu,'UnifFrictCoef = %5.3f \n',Manning_n);
    fprintf(h_mdu,'UnifFrictCoef1D = 0.023 \n');
    fprintf(h_mdu,'UnifFrictCoefLin = 0 \n');
    fprintf(h_mdu,'Umodlin = 0 \n');
    fprintf(h_mdu,'Vicouv = 0.1 \n');
    fprintf(h_mdu,'Dicouv = 0.1 \n');
    fprintf(h_mdu,'Vicoww = 1E-06 \n');
    fprintf(h_mdu,'Dicoww = 1E-06 \n');
    fprintf(h_mdu,'Vicwminb = 0 \n');
    fprintf(h_mdu,'Smagorinsky = 0.2 \n');
    fprintf(h_mdu,'Elder = 0 \n');
    fprintf(h_mdu,'Irov = 0 \n');
    fprintf(h_mdu,'wall_ks = 0 \n');
    fprintf(h_mdu,'Rhomean = 1025 \n');
    fprintf(h_mdu,'Idensform = 2 \n');
    fprintf(h_mdu,'Ag = 9.81 \n');
    fprintf(h_mdu,'TidalForcing = 1 \n');
    fprintf(h_mdu,'Doodsonstart = 55.565 \n');
    fprintf(h_mdu,'Doodsonstop = 375.575 \n');
    fprintf(h_mdu,'Doodsoneps = 0 \n');
    fprintf(h_mdu,'Salinity = 0 \n');
    fprintf(h_mdu,'InitialSalinity = 0 \n');
    fprintf(h_mdu,'Sal0abovezlev = -999 \n');
    fprintf(h_mdu,'DeltaSalinity = -999 \n');
    fprintf(h_mdu,'Backgroundsalinity = 0 \n');
    fprintf(h_mdu,'InitialTemperature = 6 \n');
    fprintf(h_mdu,'Secchidepth = 2 \n');
    fprintf(h_mdu,'Stanton = -1 \n');
    fprintf(h_mdu,'Dalton = -1 \n');
    fprintf(h_mdu,'Backgroundwatertemperature = 6 \n');
    fprintf(h_mdu,'SecondaryFlow = 0 \n');
    fprintf(h_mdu,'EffectSpiral = 0 \n');
    fprintf(h_mdu,'BetaSpiral = 0 \n');
    fprintf(h_mdu,'Temperature = 0 \n');
    fprintf(h_mdu,'VillemonteCD1 = 1 \n');
    fprintf(h_mdu,'VillemonteCD2 = 10 \n\n');

    fprintf(h_mdu,'[wind]\n');
    fprintf(h_mdu,'ICdtyp = 3 \n');
    fprintf(h_mdu,'Cdbreakpoints = %5.4f %5.4f %5.4f \n',C1,C2,C3);
    fprintf(h_mdu,'Windspeedbreakpoints = %5.1f %5.1f %5.1f \n',W1,W2,W3);
    fprintf(h_mdu,'Rhoair = 1.205 \n');
    fprintf(h_mdu,'PavBnd = %6i \n',p_n);
    fprintf(h_mdu,'PavIni = %6i \n\n',p_n);

    fprintf(h_mdu,'[waves]\n');
    fprintf(h_mdu,'Wavemodelnr = 3 \n');
    fprintf(h_mdu,'WaveNikuradse = 0.1\n');
    fprintf(h_mdu,'Rouwav = FR84 \n');
    fprintf(h_mdu,'Gammax = 1 \n\n');

    fprintf(h_mdu,'[time]\n');
    fprintf(h_mdu,'RefDate = %4i%02i%02i \n',year,month(2),day(2)');
    fprintf(h_mdu,'Tzone = 0 \n');
    fprintf(h_mdu,'DtUser = 300 \n');
    fprintf(h_mdu,'DtNodal = 21600 \n');
    fprintf(h_mdu,'DtMax = 30 \n');
    fprintf(h_mdu,'DtInit = 1 \n');
    fprintf(h_mdu,'Tunit = S \n');
    TStart = hour(2)*3600;
    fprintf(h_mdu,'TStart = %6d \n',TStart);
    TStop = time(end)-time(2)+TStart;
    fprintf(h_mdu,'TStop = %6d \n',TStop);
    fprintf(h_mdu,'StartDateTime = %4i%02i%02i%02i0000 \n',year,month(2),day(2),hour(2));
    fprintf(h_mdu,'StopDateTime = %4i%02i%02i%02i0000 \n\n',year,month(end),day(end),hour(end));

    fprintf(h_mdu,'[restart]\n');
    fprintf(h_mdu,'RestartFile = \n');
    fprintf(h_mdu,'RestartDateTime = %4i%02i%02i \n\n',year,month(2),day(2));

    fprintf(h_mdu,'[external forcing]\n');
    fprintf(h_mdu,'ExtForceFile = FlowFM.ext \n');
    fprintf(h_mdu,'ExtForceFileNew =  FlowFM_bnd.ext\n\n');

    fprintf(h_mdu,'[trachytopes]\n');
    fprintf(h_mdu,'TrtRou = N \n');
    fprintf(h_mdu,'TrtDef = \n');
    fprintf(h_mdu,'DtTrt = 60 \n\n');

    fprintf(h_mdu,'[output]\n');
    fprintf(h_mdu,'WaqInterval = 0 \n');
    fprintf(h_mdu,'Wrishp_crs = 0 \n');
    fprintf(h_mdu,'Wrishp_weir = 0 \n');
    fprintf(h_mdu,'Wrishp_gate = 0 \n');
    fprintf(h_mdu,'Wrishp_fxw = 0 \n');
    fprintf(h_mdu,'Wrishp_thd = 0 \n');
    fprintf(h_mdu,'Wrishp_obs = 0 \n');
    fprintf(h_mdu,'Wrishp_emb = 0 \n');
    fprintf(h_mdu,'Wrishp_dryarea = 0 \n');
    fprintf(h_mdu,'Wrishp_enc = 0 \n');
    fprintf(h_mdu,'Wrishp_src = 0 \n');
    fprintf(h_mdu,'Wrishp_pump = 0 \n');
    fprintf(h_mdu,'OutputDir = output \n');
    fprintf(h_mdu,'WAQOutputDir = \n');
    fprintf(h_mdu,'FlowGeomFile = \n');
    fprintf(h_mdu,'ObsFile = FlowFM_observations_obs.xyn \n');
    fprintf(h_mdu,'CrsFile = \n');
    fprintf(h_mdu,'HisFile = \n');
    fprintf(h_mdu,'HisInterval = %d \n',t_interval);
    fprintf(h_mdu,'XLSInterval = \n');
    fprintf(h_mdu,'MapFile = \n');
    fprintf(h_mdu,'MapInterval = %d \n',t_interval);
    fprintf(h_mdu,'RstInterval = 0 \n');
    fprintf(h_mdu,'S1incinterval = \n');
    fprintf(h_mdu,'MapFormat = 4 \n');
    fprintf(h_mdu,'Wrihis_balance = 1 \n');
    fprintf(h_mdu,'Wrihis_sourcesink = 0 \n');
    fprintf(h_mdu,'Wrihis_structure_gen = 1 \n');
    fprintf(h_mdu,'Wrihis_structure_dam = 1 \n');
    fprintf(h_mdu,'Wrihis_structure_pump = 0 \n');
    fprintf(h_mdu,'Wrihis_structure_gate = 0 \n');
    fprintf(h_mdu,'Wrihis_structure_weir = 1 \n');
    fprintf(h_mdu,'Wrihis_turbulence = 1 \n');
    fprintf(h_mdu,'Wrihis_wind = 1 \n');
    fprintf(h_mdu,'Wrihis_rain = 0 \n');
    fprintf(h_mdu,'Wrihis_temperature = 0 \n');
    fprintf(h_mdu,'Wrihis_heat_fluxes = 0 \n');
    fprintf(h_mdu,'Wrihis_salinity = 0 \n');
    fprintf(h_mdu,'Wrihis_density = 0 \n');
    fprintf(h_mdu,'Wrihis_waterlevel_s1 = 1 \n');
    fprintf(h_mdu,'Wrihis_waterdepth = 1 \n');
    fprintf(h_mdu,'Wrihis_velocity_vector = 1 \n');
    fprintf(h_mdu,'Wrihis_upward_velocity_component = 0 \n');
    fprintf(h_mdu,'Wrihis_sediment = 0 \n');
    fprintf(h_mdu,'Wrihis_constituents = 0 \n');
    fprintf(h_mdu,'Wrimap_waterlevel_s0 = 0 \n');
    fprintf(h_mdu,'Wrimap_waterlevel_s1 = 1 \n');
    fprintf(h_mdu,'Wrimap_velocity_component_u0 = 0 \n');
    fprintf(h_mdu,'Wrimap_velocity_component_u1 = 1 \n');
    fprintf(h_mdu,'Wrimap_velocity_vector = 1 \n');
    fprintf(h_mdu,'Wrimap_upward_velocity_component = 0 \n');
    fprintf(h_mdu,'Wrimap_density_rho = 1 \n');
    fprintf(h_mdu,'Wrimap_horizontal_viscosity_viu = 0 \n');
    fprintf(h_mdu,'Wrimap_horizontal_diffusivity_diu = 0 \n');
    fprintf(h_mdu,'Wrimap_flow_flux_q1 = 1 \n');
    fprintf(h_mdu,'Wrimap_spiral_flow = 1 \n');
    fprintf(h_mdu,'Wrimap_numlimdt = 1 \n');
    fprintf(h_mdu,'Wrimap_taucurrent = 1 \n');
    fprintf(h_mdu,'Wrimap_chezy = 0 \n');
    fprintf(h_mdu,'Wrimap_turbulence = 1 \n');
    fprintf(h_mdu,'Wrimap_wind = 1 \n');
    fprintf(h_mdu,'Wrimap_heat_fluxes = 0 \n');
    fprintf(h_mdu,'MapOutputTimeVector = \n');
    fprintf(h_mdu,'FullGridOutput = 0 \n');
    fprintf(h_mdu,'EulerVelocities = 0 \n');
    fprintf(h_mdu,'ClassMapFile = \n');
    fprintf(h_mdu,'WaterlevelClasses = 0.0 \n');
    fprintf(h_mdu,'WaterdepthClasses = 0.0 \n');
    fprintf(h_mdu,'ClassMapInterval = 0 \n');
    fprintf(h_mdu,'StatsInterval = \n');
    fprintf(h_mdu,'Writebalancefile = 0 \n');
    fprintf(h_mdu,'TimingsInterval = \n');
    fprintf(h_mdu,'Richardsononoutput = 1 \n\n');
    fclose(h_mdu);

    %% Write wave model file
    path_wave = [path_track,'\wave'];
    mkdir(path_wave);    
    copyfile([pwd,'\wave\*'],path_wave);
    copyfile([path_flow,'\',name_spider],path_wave);
    
    % Write input file (.mdw)
    name_mdw = [path_wave,'\Wave.mdw'];
    h_mdw = fopen(name_mdw,'w');
    fprintf(h_mdw,'[WaveFileInformation] \n');
    fprintf(h_mdw,'  FileVersion = 02.00 \n');

    fprintf(h_mdw,'[General]\n');
    fprintf(h_mdw,'  ReferenceDate = %4i-%02i-%02i \n',year,month(2),day(2));
    fprintf(h_mdw,'  DirConvention = nautical \n');
    fprintf(h_mdw,'  SimMode = non-stationary \n');
    fprintf(h_mdw,'  TimeStep = %d \n',t_interval/60);
    fprintf(h_mdw,'  TimeInterval = %d \n',t_interval/60);
    fprintf(h_mdw,'  OnlyInputVerify = false \n');
    fprintf(h_mdw,'  FlowBedLevel = 2 \n');
    fprintf(h_mdw,'  FlowWaterLevel = 2 \n');
    fprintf(h_mdw,'  FlowVelocity = 2 \n');
    fprintf(h_mdw,'  FlowVelocityType = depth-averaged \n');
    fprintf(h_mdw,'  FlowWind = 2 \n');
    fprintf(h_mdw,'  DirSpace = circle \n');
    fprintf(h_mdw,'  NDir = 36 \n');                             
    fprintf(h_mdw,'  StartDir = 0 \n');
    fprintf(h_mdw,'  EndDir = 360 \n');
    fprintf(h_mdw,'  NFreq = 24 \n');
    fprintf(h_mdw,'  FreqMin = 0.03 \n');
    fprintf(h_mdw,'  FreqMax = 1 \n');
    fprintf(h_mdw,'  WaterLevel = %1.4f \n',h0);
    fprintf(h_mdw,'  XVeloc = 0 \n');
    fprintf(h_mdw,'  YVeloc = 0 \n');
    fprintf(h_mdw,'  WindSpeed = 0 \n');
    fprintf(h_mdw,'  WindDir = 0 \n');
    fprintf(h_mdw,'  MeteoFile = %s \n',name_spider);

    TimeStart = hour(2)*60;
    TimeEnd = (time(end)-time(2))/60 + TimeStart;
    for t = TimeStart:t_interval/60:TimeEnd
        fprintf(h_mdw,'[TimePoint]\n');
        fprintf(h_mdw,'    Time = %e\n',t);
    end

    fprintf(h_mdw,'[Output]\n');
    fprintf(h_mdw,'  COMFile = ../dflowfm/output/FlowFM_com.nc \n');
    fprintf(h_mdw,'  WriteCOM = true \n');
    fprintf(h_mdw,'  COMWriteInterval = %d \n',t_interval/60);
    fprintf(h_mdw,'  AppendCOM = false \n');
    fprintf(h_mdw,'  MassFluxToCOM = true \n');
    fprintf(h_mdw,'  MapWriteInterval = %d \n',t_interval/60);
    fprintf(h_mdw,'  WriteTable = false \n');
    fprintf(h_mdw,'  WriteSpec1D = false \n');
    fprintf(h_mdw,'  WriteSpec2D = false \n');
    fprintf(h_mdw,'  UseHotFile = false \n');
    fprintf(h_mdw,'  TestOutputLevel = 0 \n');
    fprintf(h_mdw,'  TraceCalls = false \n');
    fprintf(h_mdw,'  MapWriteNetCDF = true \n');
    fprintf(h_mdw,'  NetCDFSinglePrecision = false \n');

    fprintf(h_mdw,'[Constants]\n');
    fprintf(h_mdw,'  WaterLevelCorrection = 0 \n');
    fprintf(h_mdw,'  Gravity = 9.81 \n');
    fprintf(h_mdw,'  WaterDensity = 1025 \n');
    fprintf(h_mdw,'  NorthDir = 90 \n');
    fprintf(h_mdw,'  MinimumDepth = 0.05 \n');

    fprintf(h_mdw,'[Processes]\n');
    fprintf(h_mdw,'  GenModePhys = 3 \n');
    fprintf(h_mdw,'  WaveSetup = false \n');
    fprintf(h_mdw,'  Breaking = true \n');
    fprintf(h_mdw,'  BreakAlpha = 1 \n');
    fprintf(h_mdw,'  BreakGamma = 0.73 \n');
    fprintf(h_mdw,'  Triads = false \n');
    fprintf(h_mdw,'  TriadsAlpha = 0.05 \n');
    fprintf(h_mdw,'  TriadsBeta = 2.5 \n');
    fprintf(h_mdw,'  BedFriction = jonswap \n');
    fprintf(h_mdw,'  BedFricCoef = 0.038 \n');
    fprintf(h_mdw,'  Diffraction = false \n');
    fprintf(h_mdw,'  DiffracSteps = 0 \n');
    fprintf(h_mdw,'  DiffracProp = true \n');
    fprintf(h_mdw,'  DiffracCoef = 0.2 \n');
    fprintf(h_mdw,'  WindGrowth = true \n');
    fprintf(h_mdw,'  Quadruplets = true \n');
    fprintf(h_mdw,'  WhiteCapping = Komen \n');
    fprintf(h_mdw,'  Refraction = true \n');
    fprintf(h_mdw,'  FreqShift = true \n');
    fprintf(h_mdw,'  WaveForces = dissipation 3d \n');

    fprintf(h_mdw,'[Numerics]\n');
    fprintf(h_mdw,'  DirSpaceCDD = 0.5 \n');
    fprintf(h_mdw,'  FreqSpaceCSS = 0.5 \n');
    fprintf(h_mdw,'  RChHsTm01 = 0.02 \n');
    fprintf(h_mdw,'  RChMeanHs = 0.02 \n');
    fprintf(h_mdw,'  RChMeanTm01 = 0.02 \n');
    fprintf(h_mdw,'  PercWet = 98 \n');
    fprintf(h_mdw,'  MaxIter = 15 \n');

    fprintf(h_mdw,'[Boundary]\n');
    fprintf(h_mdw,'  Name = Boundary \n');
    fprintf(h_mdw,'  Definition = xy-coordinates \n');
    fprintf(h_mdw,'  StartCoordX = -98.3189698 \n');
    fprintf(h_mdw,'  EndCoordX = -86.4447522 \n');
    fprintf(h_mdw,'  StartCoordY = 22.3209980 \n');
    fprintf(h_mdw,'  EndCoordY = 30.9557683 \n');
    fprintf(h_mdw,'  SpectrumSpec = parametric \n');
    fprintf(h_mdw,'  SpShapeType = Jonswap \n');
    fprintf(h_mdw,'  PeriodType = peak \n');
    fprintf(h_mdw,'  DirSpreadType = power \n');
    fprintf(h_mdw,'  PeakEnhanceFac = 3.3000000e+000 \n');
    fprintf(h_mdw,'  WaveHeight = 0.0000000e+000 \n');
    fprintf(h_mdw,'  Period = 1.0000000e+000 \n');
    fprintf(h_mdw,'  Direction = 0.0000000e+000 \n');
    fprintf(h_mdw,'  DirSpreading = 4.0000000e+000 \n');

    fprintf(h_mdw,'[Domain]\n');
    fprintf(h_mdw,'  Grid = Grid_entire_net.grd \n');
    fprintf(h_mdw,'  BedLevel = Bathymetry_entire.dep \n');
    fprintf(h_mdw,'  Output = true \n');
    fprintf(h_mdw,'[Domain]\n');
    fprintf(h_mdw,'  Grid = Grid_nesting_net.grd \n');
    fprintf(h_mdw,'  BedLevel = Bathymetry_nesting.dep \n');
    fprintf(h_mdw,'  NestedInDomain = 1 \n');
    fprintf(h_mdw,'  Output = true \n');
    fprintf(h_mdw,'[Domain]\n');
    fprintf(h_mdw,'  Grid = Grid_nesting01.grd \n');
    fprintf(h_mdw,'  BedLevel = Bathymetry_nesting01.dep \n');
    fprintf(h_mdw,'  NestedInDomain = 2 \n');
    fprintf(h_mdw,'  Output = true \n');
    fprintf(h_mdw,'[Domain]\n');
    fprintf(h_mdw,'  Grid = Grid_nesting02.grd \n');
    fprintf(h_mdw,'  BedLevel = Bathymetry_nesting02.dep \n');
    fprintf(h_mdw,'  NestedInDomain = 2 \n');
    fprintf(h_mdw,'  Output = true \n');
    fprintf(h_mdw,'[Domain]\n');
    fprintf(h_mdw,'  Grid = Grid_nesting03.grd \n');
    fprintf(h_mdw,'  BedLevel = Bathymetry_nesting03.dep \n');
    fprintf(h_mdw,'  NestedInDomain = 2 \n');
    fprintf(h_mdw,'  Output = true \n');
    fprintf(h_mdw,'[Domain]\n');
    fprintf(h_mdw,'  Grid = Grid_nesting04.grd \n');
    fprintf(h_mdw,'  BedLevel = Bathymetry_nesting04.dep \n');
    fprintf(h_mdw,'  NestedInDomain = 2 \n');
    fprintf(h_mdw,'  Output = true \n');
    fclose(h_mdw);

    %% Write DIMR file
    name_dimr = [path_track,'\dimr_config.xml'];
    h_dimr = fopen(name_dimr ,'w');
    fprintf(h_dimr,'<?xml version="1.0" encoding="utf-8" standalone="yes"?> \n');
    fprintf(h_dimr,'<dimrConfig xmlns="http://schemas.deltares.nl/dimr" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://schemas.deltares.nl/dimr http://content.oss.deltares.nl/schemas/dimr-1.2.xsd"> \n');
    fprintf(h_dimr,'  <documentation>\n');
    fprintf(h_dimr,'    <fileVersion>1.2</fileVersion> \n');
    fprintf(h_dimr,'    <createdBy>Deltares, Coupling Team</createdBy> \n');
    fprintf(h_dimr,'    <creationDate>2024-01-10T21:20:59.5818399Z</creationDate> \n');
    fprintf(h_dimr,'  </documentation>\n');

    fprintf(h_dimr,'  <control> \n');
    fprintf(h_dimr,'    <parallel> \n');
    fprintf(h_dimr,'      <startGroup> \n');
    fprintf(h_dimr,'        <time>0 %d %7d</time> \n',t_interval,TStop);
    fprintf(h_dimr,'        <start name="Wave" /> \n');
    fprintf(h_dimr,'      </startGroup> \n');
    fprintf(h_dimr,'      <start name="FlowFM" /> \n'); 
    fprintf(h_dimr,'    </parallel> \n');
    fprintf(h_dimr,'  </control> \n');

    fprintf(h_dimr,'  <component name="Wave"> \n');
    fprintf(h_dimr,'    <library>wave</library> \n');
    fprintf(h_dimr,'    <workingDir>wave</workingDir> \n');
    fprintf(h_dimr,'    <inputFile>Wave.mdw</inputFile> \n');
    fprintf(h_dimr,'  </component> \n');

    fprintf(h_dimr,'  <component name="FlowFM"> \n');
    fprintf(h_dimr,'    <library>dflowfm</library> \n');
    fprintf(h_dimr,'    <workingDir>dflowfm</workingDir> \n');
    fprintf(h_dimr,'    <inputFile>FlowFM.mdu</inputFile> \n');
    fprintf(h_dimr,'  </component> \n');
    fprintf(h_dimr,'</dimrConfig> \n');
    fclose(h_dimr);

end

