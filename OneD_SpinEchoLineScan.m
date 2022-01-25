%% Spin Echo Line Scan
% collects just the desired read line to speed up functionality

clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-5; 
gamma = 2*pi*42.577*10^6;

%load the voxel model
load('PD.mat');
load('T1.mat');
load('T2.mat');

SCALE_FACTOR_x = 16;
SCALE_FACTOR_y = 1;

%Allocate the memory needed
nTimeSteps  = 700*(10^-5/dt);%*48;
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(2,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points

xStepsSim  = size(T1,1);   %Number of simulated "spins" in the x directions 
yStepsSim  = size(T1,2);   %Number of simulated "spins" in the y directions 
zStepsSim  = 1;            %Number of simulated "spins" in the z directions 

dX = (4.0e-3)/SCALE_FACTOR_x;            %Distance between simulated "spins" in the x directions  [meter]
dY = (4.0e-3)/SCALE_FACTOR_y;            %Distance between simulated "spins" in the y directions  [meter]
dZ = (1.0e-4);            %Distance between simulated "spins" in the z directions  [meter]

xStepsImg = 48;
yStepsImg = 48;
zStepsImg = 1;

% %% SELECT A LINE %% %
% LINE_SELECTED = 25;

% Ask user for a number of steps to take.
% defaultValues = {'25'};
% titleBar = 'What line would you like to select';
% userPrompt = {sprintf('Select a line from %d-%d: ',0,xStepsImg)};
% caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
% if isempty(caUserInput),return,end; % stop if clicked cancel
% 
% LINE_SELECTED = round(str2num(caUserInput{1}));

LINE_SELECTED = 25;

%check what line based on centerline
if (mod(xStepsImg,2) == 0) %even
    center_line = (xStepsImg/2)+1;
    
else %odd
    center_line = ((xStepsImg-1)/2)+1;
end

actual_line = LINE_SELECTED - center_line; %need to figure out if it will mov up or down if neg or positive

%Find paramters for gradients 
TE = 6*1.0e-3; %3ms for TE/2
data_points = xStepsImg;
FOV_r=192*(10^-3);
delta_r = FOV_r/data_points ;
dwellt = 10*1.0e-6;
t_acq = dwellt*data_points;
G_r = ((2*pi)/gamma) * (1/(t_acq*delta_r)); %(2*pi)/gamma)*(1/(delta_r*t_acq))
BW_p = (6/0.001); %2n/tau, manually just made 3 for example
G_s = BW_p*(2*pi/(gamma*0.005));%slice thickness is 5mm 


% 3D positions in space
pos = zeros(3,xStepsSim*SCALE_FACTOR_x,yStepsSim*SCALE_FACTOR_y,zStepsSim);
for k=1:xStepsSim*SCALE_FACTOR_x
    for j=1:yStepsSim*SCALE_FACTOR_y
        for i=1:zStepsSim
            pos(1,k,j,i) = (k-xStepsSim*SCALE_FACTOR_x/2)*dX;
            pos(2,k,j,i) = (j-yStepsSim*SCALE_FACTOR_y/2)*dY;
            pos(3,k,j,i) = (i-zStepsSim/2)*dZ;
        end
    end
end



%Generates the time line for sequence plotting
for i=1:nTimeSteps 
    time(i)    = i*dt;                       %Time in seconds
end

% %% RF %% %
% RF Excitation Waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = apodize_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

%RF Refocusing pules
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
startR = (TE/2);

rfPulseR = apodize_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
%Find the line that was selected and edit phase ramp/shift
SELECTED_LINE = exp(1j*BW_p*(actual_line)*(1:(length(rfPulseR))));

%line scan selection
rfPulseR = rfPulseR.*SELECTED_LINE;
rfPulse(round((startR*(1/dt)))+(rfStepsR)-1) = rfPulseR;

% %% Gradients  %% %
% Generate read and Phase encoding for signal line
% Generate the phase encoding gradient
timestepsPHASEPrephaseEvent = round(((TE/dt)-(74*10^-5)/dt):((TE/dt)-(25*10^-5)/dt));
gradAmp(2,timestepsPHASEPrephaseEvent) = -G_r * (48/50) /2;  %Y gradients in Tesla per meter

% Read + ADC
timestepsREADOUTPrephaseEvent = round(((TE/dt)-(74*10^-5)/dt):((TE/dt)-(25*10^-5)/dt));

% timestepsADCEvent = 576:(576+(xStepsImg)-1);
% timestepsREADOUTEvent = 576:(576+(xStepsImg)-1);


timestepsADCEvent = round((((TE/dt)-(25*10^-5)/dt)+(10^-5/dt)):((((TE/dt)-(25*10^-5)/dt)+(10^-5/dt))+(480*10^-6)/dt)-1);
timestepsREADOUTEvent = round((((TE/dt)-(25*10^-5)/dt)+(10^-5/dt)):((((TE/dt)-(25*10^-5)/dt)+(10^-5/dt))+(480*10^-6)/dt)-1);

gradAmp(1,timestepsREADOUTPrephaseEvent) =  (-G_r * (48/50) / 48); % X read prephase
%change this for only line increments, for test lets do 1 point
%should prephase the first point to b isocenter

readout_steps = 1:(xStepsImg);
readout_steps=readout_steps(repmat(1:size(readout_steps,1),1,1),repmat(1:size(readout_steps,2),(10^-5/dt),1))

% adc(1,timestepsREADOUTEvent) = 1:48; % adc event
gradAmp(1,timestepsREADOUTEvent) = G_r/48; %X read out event

% adc(1,576+LINE_SELECTED-1) = LINE_SELECTED; % adc event

adc(1,timestepsADCEvent) = 1:(xStepsImg); % adc event




%%% ADDDING SPOILING GRADIENTS
% n=3;
% G_spoil = (n*((2*pi)/(gamma*0.005)))/(dt*50);
% gradAmp(1,(576+(xStepsImg)):((576+(xStepsImg))+49)) =  -G_r * (48/50) / 2;


% %% Spatial Encoding

gradAmp(3,round(1:(pulsedurE/(dt)))) =  G_s; %Z gradients in Tesla per meter
gradAmp(3,round((pulsedurE/(dt))+1):((pulsedurR/(dt))+((pulsedurR/(dt))/2))) =  -G_s; %Z gradients in Tesla per meter

gradAmp(1,round((startR*(1/dt)))+(rfStepsR)) =  G_s; %Z gradients in Tesla per meter orig
gradAmp(1,round((startR*(1/dt))-((pulsedurR/(dt))/2)+(1:((pulsedurR/(dt))/2)))) =  3.2*G_s; %Z gradients in Tesla per meter
gradAmp(1,round((startR*(1/dt)))+(pulsedurR/(dt)+(1:((pulsedurR/(dt))/2)))) =  3.2*G_s; %Z gradients in Tesla per meter


% %% PLOTTING %%
figure
subplot(6,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('Full Sequence: RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
subplot(6,1,2); plot(time,gradAmp(1,:),'r-','LineWidth',2);title('Full Sequence: Read Gradient'); ylim([-0.06 0.05]);
xlabel('time (s)'), ylabel('G_{r}(T/m)');grid on;
subplot(6,1,3); plot(time,gradAmp(2,:),'g-','LineWidth',2);title('Full Sequence: Phase Gradient');ylim([-0.06 0.05]);
xlabel('time (s)'), ylabel('G_{p}(T/m)');grid on;
subplot(6,1,4); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('Full Sequence: Slice Select Gradient');ylim([-0.06 0.05])
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

subplot(6,1,5); plot(time,adc(1,:),'k-','LineWidth',2);title('ADC_{x} Events');ylim([-20 60])
xlabel('time (s)'), ylabel('ADC Event Count');grid on;
subplot(6,1,6); plot(time,adc(2,:),'k-','LineWidth',2);title('ADC_{y} Event');ylim([-20 60])
xlabel('time (s)'), ylabel('ADC Event Count');grid on;

%% Perform seuquence
kSpace  = zeros(yStepsImg,yStepsImg); %vector since 1d line reading
durPE = 0.5e-3;
gradPEmax = pi * (yStepsSim-1)/(gamma*FOV_r) * 1/durPE;


T1_new = T1(repmat(1:size(T1,1),SCALE_FACTOR_x,1),repmat(1:size(T2,2),SCALE_FACTOR_y,1));
T2_new = T2(repmat(1:size(T2,1),SCALE_FACTOR_x,1),repmat(1:size(T2,2),SCALE_FACTOR_y,1));

tic

%Cycle through all kspace for each row
for trIndex=1:yStepsImg
    disp(trIndex);  
    
%     update the phase endcoding gradient
    gradAmp(2,timestepsREADOUTPrephaseEvent) =  (1-(trIndex-1)/((yStepsSim-1)/2)) * gradPEmax; %Y gradients in Tesla per meter
    
    % Set up spin selection vectors
    
 x_Steps = xStepsSim*SCALE_FACTOR_x;
 y_Steps = yStepsSim*SCALE_FACTOR_y;
    
k=1:x_Steps;  
j=1:y_Steps;
counter = 1;
    for i=1:zStepsSim

        % Set up inital location
        r = reshape(pos(:,k,j,i),[3,x_Steps*y_Steps]);
        mT = zeros(x_Steps*y_Steps,1);
        mZ = ones(x_Steps*y_Steps,1);

        % Find inital magnetization
        dB0 = ((gradAmp(:,1)'*r)')+2.5e-7*(rand(x_Steps*y_Steps,1)-1);%+2.5e-7*(rand(1,nSpins)*2-1); 
        [mT,mZ] =  bloch(dt, dB0,rfPulse(1), reshape(T1_new(k,j),[x_Steps*y_Steps,1]),...
                                             reshape(T2_new(k,j),[x_Steps*y_Steps,1]),...
                                             mT, mZ);  

        % Find magnetization for rest of sequence                              
        for t=2:nTimeSteps 
            dB0 = ((gradAmp(:,t)'*r)')+2.5e-7*(rand(x_Steps*y_Steps,1)-1);
            [mT,mZ] =  bloch(dt, dB0,rfPulse(t), reshape(T1_new(k,j),[x_Steps*y_Steps,1]),...
                                                 reshape(T2_new(k,j),[x_Steps*y_Steps,1]), ...
                                                 mT, mZ);

            % Encode data for ADC readout measurements
            if(adc(1,t)>0)
                PDt_new = PD(repmat(1:size(PD,1),SCALE_FACTOR_x,1),repmat(1:size(PD,2),SCALE_FACTOR_y,1));
                PDt = PDt_new(k,j);
                
                kSpace(round(adc(1,t)),trIndex) = kSpace(round(adc(1,t)),trIndex)...
                                                  + sum(mT.*reshape(PDt(k,j),[x_Steps*y_Steps,1]));

            end %end of if
        end  %end of time loop                     
    end %end z steps
end %end of TR loop
t = toc;

%% Reconstructing the image

figure;
Ispace = fftshift(ifft2(fftshift(kSpace)));
subplot(2,1,1);imagesc(abs(Ispace));
title('Magnitude Image for Center most Line (line 25)');
subplot(2,1,2);imagesc(log(abs((kSpace))),'CDataMapping','scaled'); 
title('Absolute kSpace');

save('LineImage.mat','Ispace')