%% Spin Echo Line Scan
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-5; 
gamma = 2*pi*42.577*10^6;

%load the voxel model
load('PD.mat');
load('T1.mat');
load('T2.mat');

SCALE_FACTOR_x = 8;
SCALE_FACTOR_y = 1;

%Allocate the memory needed
nTimeSteps  = 700;%*48;
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


% Generate and Display MRI Sequence
TE = 6*1.0e-3; %3ms for TE/2
data_points = xStepsImg;
FOV_r=192*(10^-3);
delta_r = FOV_r/data_points ;
dwellt = 10*1.0e-6;
t_acq = dwellt*data_points;
G_r = ((2*pi)/gamma) * (1/(t_acq*delta_r)); %(2*pi)/gamma)*(1/(delta_r*t_acq))
BW_p = (6/0.001); %2n/tau
G_s = BW_p*(2*pi/(gamma*0.005));%slice thickness is 5mm 

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
%line scan selection
rfPulseR = rfPulseR.*exp(1j*24000*(1:100));
rfPulse(round((startR*10^5))+(rfStepsR)-1) = rfPulseR;

% %% Gradients  %% %

% Generate the phase encoding gradient
timestepsPHASEPrephaseEvent = round(((TE/dt)-(74*10^-5)/dt):((TE/dt)-(25*10^-5)/dt));
gradAmp(2,timestepsPHASEPrephaseEvent) = G_r * (48/50) /2;  %Y gradients in Tesla per meter

% Read + ADC
timestepsREADOUTPrephaseEvent = round(((TE/dt)-(74*10^-5)/dt):((TE/dt)-(25*10^-5)/dt));
timestepsADCEvent = 576:(576+(xStepsImg)-1);
timestepsREADOUTEvent = 576:(576+(xStepsImg)-1);

gradAmp(1,timestepsREADOUTPrephaseEvent) =  -G_r * (48/50) / 2; % X read prephase
adc(1,timestepsADCEvent) = 1:(xStepsImg); % adc event
gradAmp(1,timestepsREADOUTEvent) = G_r; %X read out event
adc(2,timestepsADCEvent) = ones(1,xStepsImg);

%%% ADDDING SPOILING GRADIENTS
% n=3;
% G_spoil = (n*((2*pi)/(gamma*0.005)))/(dt*50);
% gradAmp(1,round((startR*10^5)-((pulsedurR/(dt))/2)+(1:((pulsedurR/(dt))/2)))) =  G_spoil;
% 
gradAmp(1,(576+(xStepsImg)):((576+(xStepsImg))+49)) =  -G_r * (48/50) / 2;
% gradAmp(2,(576+(xStepsImg)):((576+(xStepsImg))+49)) =  -G_r * (48/50) / 2;
% %%%

% %% Spatial Encoding

gradAmp(3,round(1:(pulsedurE/(dt)))) =  G_s; %Z gradients in Tesla per meter
gradAmp(3,round((pulsedurE/(dt))+1):((pulsedurR/(dt))+((pulsedurR/(dt))/2))) =  -G_s; %Z gradients in Tesla per meter

gradAmp(1,round((startR*10^5))+(rfStepsR)) =  G_s; %Z gradients in Tesla per meter orig
gradAmp(1,round((startR*10^5)-((pulsedurR/(dt))/2)+(1:((pulsedurR/(dt))/2)))) =  2*G_s; %Z gradients in Tesla per meter
gradAmp(1,round((startR*10^5))+(pulsedurR/(dt)+(1:((pulsedurR/(dt))/2)))) =  2*G_s; %Z gradients in Tesla per meter



%% Create entire TR Sequence
% TOTAL_REPETITIONS = 48;
% gradPEmax = pi * (yStepsSim-1)/(gamma*FOV_r) * 1/0.5e-3;
% for i = 1:47
%     nStart =(700*i)+1;
%     nEnd   =(700*i)+700;
%     
%     rfPulse(1,nStart:nEnd) = rfPulse(1,1:700);
% 
%     gradAmp(1,nStart:nEnd) = gradAmp(1,1:700);
%     gradAmp(2,nStart-1+timestepsPHASEPrephaseEvent) = (1-(i+1-1)/((yStepsSim-1)/2)) * gradPEmax; %Y gradients in Tesla per meter
%     gradAmp(2,nStart-1+((576+(xStepsImg)):((576+(xStepsImg))+49))) = -1* (1-(i+1-1)/((yStepsSim-1)/2)) * gradPEmax; %Y gradients in Tesla per meter
% 
%     gradAmp(3,nStart:nEnd) = gradAmp(3,1:700); 
%     
%     adc(1,nStart:nEnd)     = adc(1,1:700);
%     adc(2,nStart:nEnd)     = adc(2,1:700)*(i+1);
% 
% end
% 
% %plot the compleate sequence
% figure
% subplot(6,1,1); plot(time(1:700),rfPulse(1:700),'k-','LineWidth',2);title('TR Cycle: RF Pulse'); 
% xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;
% subplot(6,1,2); plot(time(1:700),gradAmp(1,(1:700)),'r-','LineWidth',2);title('TR Cycle: Read Gradient'); ylim([-0.06 0.05]);
% xlabel('time (s)'), ylabel('G_{r}(T/m)');grid on;
% subplot(6,1,3); plot(time(1:700),gradAmp(2,(1:700)),'g-','LineWidth',2);title('TR Cycle: Phase Gradient');ylim([-0.06 0.05]);
% xlabel('time (s)'), ylabel('G_{p}(T/m)');grid on;
% subplot(6,1,4); plot(time(1:700),gradAmp(3,(1:700)),'b-','LineWidth',2);title('TR Cycle: Slice Select Gradient');ylim([-0.06 0.05])
% xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;
% 
% subplot(6,1,5); plot(time(1:700),adc(1,(1:700)),'k-','LineWidth',2);title('ADC_{x} Event');ylim([-20 60])
% xlabel('time (s)'), ylabel('ADC Event Count');grid on;
% subplot(6,1,6); plot(time(1:700),adc(2,(1:700)),'k-','LineWidth',2);title('ADC_{x} Event');ylim([-20 60])
% xlabel('time (s)'), ylabel('ADC Event Count');grid on;




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
kSpace  = zeros(xStepsImg,yStepsImg);
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
%     gradAmp(2,(576+(xStepsImg)):((576+(xStepsImg))+50)) =   (1-(trIndex-1)/((yStepsSim-1)/2)) * gradPEmax;
    
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

%                 kSpace(round(adc(1,t)),round(adc(2,t))) = kSpace(round(adc(1,t)),round(adc(2,t)))...
%                                                 +sum(mT.*reshape(PDt(k,j),[x_Steps*y_Steps,1]));
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
