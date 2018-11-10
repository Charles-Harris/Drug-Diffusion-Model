function hw6

clear all;
clc;

% Initializing Variables

A = 0.25; % Maxium alpha for stability
dH = 0.5; % Delta h in mm
D = 6.8; % Diffusion constant in mm^2/hr
C = 2.3; % Concentration in ug/g

% Calculating delta T
dT = (A*(dH)^2)/D;

% Initializing 60 mm x 60 mm grid

% Sets the bounds of the grid
negativeXLimit = -30;
positiveXLimit = 30;
negativeYLimit = -30;
positiveYLimit = 30;

% Since dX = dY = dH
xRange = negativeXLimit:dH:positiveXLimit;
yRange = negativeYLimit:dH:positiveYLimit;

% Time variables
time = 0:dT:120; %dT is change in time

% This determines the index of the desired time values
concentrationTime15 = find(abs(time-15)<(dT/2));
concentrationTime30 = find(abs(time-30)<(dT/2));

% Implant placed at origin -> concentration at origin must = C

% Find Origin Location

% Add 1 because length is 121 and an integer result is desired
xOrigin = (1+length(xRange))/2;
yOrigin = (1+length(yRange))/2;

% Create a storage unit for concentration values on grid
concentrationStorage = zeros(length(xRange),length(yRange));

% Declaring variables for returning storage values at key points

% Origin
concentrationAtOriginX = find(xRange==0);
concentrationAtOriginY = find(yRange == 0);

% 6mm Away
storage6 = zeros(length(time),1);
concentrationAt6 = find(xRange==6);

% 12mm Away
storage12 = zeros(length(time),1);
concentrationAt12 = find(xRange==12);

% Populate storage unit with values

% Everything thats not the origin or the edges

% Indicies are established so as not to affect constant regions
i = 2:length(xRange)-1;
j = 2:length(yRange)-1;

for timeValue = 2:length(time)
    
    concentrationStorage(i,j) = concentrationStorage(i,j)*(1-4*A)+A*(concentrationStorage(i+1,j)+concentrationStorage(i-1,j)+concentrationStorage(i,j+1)+concentrationStorage(i,j-1));

    % Origin
    concentrationStorage(xOrigin,yOrigin) = C;

    % Edges (simulate infinite distance)

    % Bottom edge
    concentrationStorage(:,1) = 0;
    
    % Left Edge
    concentrationStorage(1,:)= 0;
    
    %Right Edge
    concentrationStorage(121,:)=0;
    
    %Top Edge
    concentrationStorage(:,121) = 0;
    
    % Saves concentration at grid points t =15 and t=30
    if timeValue == concentrationTime15
        concentrationStorage15 = concentrationStorage;
    end
    
    if timeValue == concentrationTime30
        concentrationStorage30 = concentrationStorage;
    end
    
    % Saves concentrations at 6 and 12 mm away from the implant
    storage6(timeValue) = concentrationStorage(concentrationAt6, concentrationAtOriginY);
    storage12(timeValue) = concentrationStorage(concentrationAt12,concentrationAtOriginY);
    
end
disp(storage6((length(time)+1)/2))
time120 = length(time);
% Outputs concentrations at t =120 for desired distances
disp('The Concentration 6 mm away at time = 120 hours is:')
disp(storage6(time120))
disp('The Concentration 12 mm away at time = 120 hours is:')
disp(storage12(time120))

% Graphing

% Cross sections

figure;

hold on;
% concentrationAtOriginY effectively sets y=0
% ':' covers the whole x range from -30 to 30
crossOne = plot(xRange,concentrationStorage15(:,concentrationAtOriginY),'r');
crossTwo = plot(xRange,concentrationStorage30(:,concentrationAtOriginY),'y');
crossThree = plot (xRange,concentrationStorage(:,concentrationAtOriginY),'k');
hold off;
ylabel('Concentration (micrograms)');
title('Concentration of Drug (micrograms) vs Distance (millimeters) Cross Section');
xlabel('Distance (millimeters)');
legend([crossOne crossTwo crossThree], {'t = 15','t = 30', 't=120'});


% Concentration Ratio vs Time
figure;
ratioPlot = plot(time, storage12./storage6,'k');
ylabel('Concentration Ratio (Dimensionless)');
xlabel('Time (hours)');
title('Concentration Ratio(dimensionless) vs Time (hours)');
legend([ratioPlot],{'Ratio of C(12mm) to C(6mm)'});

% Individual Concentrations vs Time
figure;
hold on;
c12 = plot(time,storage12, 'r');
c6 = plot(time,storage6,'y');
hold off;
ylabel('Concentration (micrograms)');
xlabel('Time(hours)');
title('Concentration (micrograms) vs Time (Hours)');
legend([c6 c12], {'C(6mm)', 'C(12mm)'});

% Mesh
figure;
mesh(xRange,yRange,concentrationStorage);
ylabel('X Distance From Origin (mm)');
xlabel('Y Distance From Origin(mm)');
zlabel('Concentration (micrograms)');
title('Concentration (micrograms) vs Distance (mm) at time = 120 hours');
end