function out = preCalcVisc(airfoil,chord)

NACA = split(airfoil,"_");
airfoilCoord = createNACA4(NACA{2},99);
foilcords = [airfoilCoord(1:100,1), airfoilCoord(1:100,2), airfoilCoord(100:199,1), airfoilCoord(100:199,2)];

%% xfoil polars for different reynolds
alpha = -10:1:20;

V = [30 80];
rho = 1.225; % kg/m³
mu = 1.81*10^-5; % kg/(m·s)

ReA = rho*V.*chord/mu;

for i=1:length(ReA)
    Mach = V(i)/340.29;
    Re = ReA(i);

    xf = XFOIL;
    xf.KeepFiles = false; % Set it to true to keep all intermediate files created (Airfoil, Polars, ...)
    xf.Visible = false;    % Set it to false to hide XFOIL plotting window
    xf.ID = V(i);
    xf.Actions = {};
    xf.Airfoil =  Airfoil.useCord(foilcords); % Define Airfoil
    %Add five filtering steps to smooth the airfoil coordinates and help convergence
    xf.addFiltering(5);
    %Switch to OPER mode, and set Reynolds = 3E7, Mach = 0.1
    xf.addOperation(Re,Mach);
    %Set maximum number of iterations
    xf.addIter(200)
    %Initializate the calculations
    xf.addAlpha(0,true);
    %Create a new polar
    xf.addPolarFile(sprintf('Polar_%d.txt', V(i)));
    % xf.addPolarFile(sprintf('Polar.txt'));
    %Calculate a sequence of angle of attack, from 0 to 25 degrees, step size of 0.1 degrees
    xf.addAlpha(alpha);
    %Close the polar file
    xf.addClosePolarFile;
    %And finally add the action to quit XFOIL
    xf.addQuit;
    % run XFOIL
    xf.run
    % Wait up to X seconds for it to finish... 
    finished = xf.wait(20);
    xf.kill;
    xf.readPolars;
    if numel(xf.Polars{1}.Alpha) > 1
        out.cl(i,:) = interp1(xf.Polars{1}.Alpha,xf.Polars{1}.CL,alpha,'linear');
        out.cd(i,:) = interp1(xf.Polars{1}.Alpha,xf.Polars{1}.CD,alpha,'linear');

        out.alpha(i,:) = alpha;
        out.ReINT(i,:) = Re*ones(1,length(alpha));
        out.Re(i) = Re;
        out.maxCl(i) = max(xf.Polars{1}.CL);
        
        out.cl(i,isnan(out.cd(i,:))) = 0;
        out.cd(i,isnan(out.cd(i,:))) = max(out.cd(i,:));
    end
end

out.chordMean = chord;




