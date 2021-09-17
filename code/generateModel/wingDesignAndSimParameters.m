function [wingDesign, simParam] = wingDesignAndSimParameters(morphingWing, airfoil, plt)


%% define wing design

% general wing parameters
wingDesign.airfoil = airfoil; % choose airfoil: coordinates are loaded from file in folder 'airfoils' 
wingDesign.chord = 0.5; % chord length
wingDesign.span = 10; % span 

% spars, ribs
wingDesign.nRibs = 9;               % number of rib pairs
wingDesign.ribActiveWidth = 0.05;   % rib pair width (width between 2 compliant ribs of a rib pair)
wingDesign.ribWidth = 0.01;         % single rib width
wingDesign.rearSpar = 0.3;          % rear spar X position, relative to chord length, measured from leading edge
wingDesign.frontSparX = 2/3;        % front spar X position, relative to distance between leading edge and rear spar. measured from leading edge 

% layup thickness
wingDesign.t_spar = 0.002;      % thickness front and rear spar at wing root
wingDesign.t_skin = 0.001;      % thickness  
wingDesign.flexSkinT = 0.0006;  % thickness of flexible skin bottom wing surface
wingDesign.rigidRibT = 0.001;   % thickness of rigid ribs
wingDesign.actuationT = 0.002;  % 0.001; % thickness of actuation mounting

% material properties
wingDesign.matProp.E = [70e9, 70e9, 1e7];       % young's modulus for 3 materials: wing skin, internal morphing structure, flexible skin bottom wing surface
wingDesign.matProp.rho = [1600, 1600, 1600];    % density
wingDesign.matProp.nu = [0.4, 0.4, 0.4];        % Poisson's ratio

% morphing parameters
wingDesign.nRibsC = 5;              % number of compliant morphing rib pairs
wingDesign.corrugationStart = 0.4;  % corrugation (flexible skin on lower surface) start X position, relative to distance between rear spar and trailing edge, measured from rear spar
wingDesign.corrugationEnd = 0.5;    % corrugation (flexible skin on lower surface) end X position, relative to distance between rear spar and trailing edge, measured from rear spar
wingDesign.actuatorFrontX = 0.05;   % actuator mounting front (at rear spar) X position, relative to distance between rear spar and trailing edge, measured from rear spar
wingDesign.actuatorFrontY = 0.2;    % actuator mounting front (at rear spar) Y position, relative to airfoil thickness, measured from lower side
wingDesign.actuatorRearX = 0.7;     % actuator mounting rear X position, relative to distance between rear spar and trailing edge, measured from rear spar
wingDesign.actuatorRearY = 0.3;     % actuator mounting rear Y position, relative to airfoil thickness, measured from lower side
% voronoi parameters: Voronoi X (1:12), Voronoi Y (13:24), Voronoi thickness (25:36)
wingDesign.voronoiXY = [1.10000000000000;0.391969163035072;0.179488214230541;0.160414284166679;0.0841243668778740;0.374190194558112;-0.0244206167537660;1.00479813400891;0.917933187941572;0.475269448973662;0.792382179944412;-0.0423462835530037;0.0684637720063988;0.169399481966369;0.0910962779278634;0.184037684888475;0.110522344877812;0.160485547131127;0.226772667987358;0;0.117090964059230;0.157742185996391;0.246589004535733;0];
wingDesign.voronoiT = 0.0004;       % thickness morphing internal rib structure

% plotting options
wingDesign.plot = plt; % plot FE mesh


%% define main FSI simulation parameters

% simParam.alpha = 6; % base angle of attack
% simParam.Velocity = 50; % 
simParam.actForce = 1000; % actuator force
simParam.nmodes = 8; % number of structure vibration modes. In case of non moprhing wing, 2 additional vibration modes are added
simParam.addMM = 1; % add morphing modes to vibration modes -> method to reduce number of vibration modes: add moprhing deformation "modes" to vibration modes

% structural mesh
simParam.elem_max_sizeC = 0.01; % maximum element size FE mesh chordwise direction
simParam.elem_max_sizeC_LE = 0.01; % maximum element size FE mesh chordwise direction around leading edge
simParam.elem_max_sizeS = 0.16; % maximum element size FE mesh spanwise direction

% aero mesh
simParam.num_airfoil_nodes_panel = 20; % aero mesh number of panels chordwise direction
simParam.n_seg_PM = 8; % aero mesh number of panels spanwise direction
simParam.n_seg_LLT = 9; % number of panels lifting line for drag calcualtion
simParam.rho = 1.225; % rho air
simParam.simtol = 10^-3;%10^-4; % FSI convergence criteria on relative pressure coefficient and wing deformation


%% non-morphing wing: overwrite these parameters to generate a non-morphing wing
if ~morphingWing
    wingDesign.nRibsC = 0; % number of compliant rib pairs
    wingDesign.matProp.E(3) = wingDesign.matProp.E(1); 
    simParam.addMM = 0; % add morphing modes to vibration modes
end
