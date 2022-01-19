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


%Allocate the memory needed
nTimeSteps  = 700;
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points


xSteps  = size(T1,1);   %Number of simulated "spins" in the x directions 
ySteps  = size(T1,2);   %Number of simulated "spins" in the y directions 
zSteps  = 1;            %Number of simulated "spins" in the z directions 

dX = 4.0e-3;            %Distance between simulated "spins" in the x directions  [meter]
dY = 4.0e-3;            %Distance between simulated "spins" in the y directions  [meter]
dZ = 1.0e-4;            %Distance between simulated "spins" in the z directions  [meter]


% 3D positions in space
pos = zeros(3,xSteps,ySteps,zSteps);
for k=1:xSteps
    for j=1:ySteps
        for i=1:zSteps
            pos(1,k,j,i) = (k-xSteps/2)*dX;
            pos(2,k,j,i) = (j-ySteps/2)*dY;
            pos(3,k,j,i) = (i-zSteps/2)*dZ;
        end
    end
end

%Generates the time line for sequence plotting
for i=1:nTimeSteps 
    time(i)    = i*dt;                       %Time in seconds
end


% Generate and Display MRI Sequence

TE = 6*1.0e-3; %3ms
data_points = 48;
FOV_r=192*(10^-3);
delta_r = FOV_r/data_points ;
dwellt = 10*10^-6;
t_acq = dwellt*data_points;
G_r = ((2*pi)/gamma) * (1/(t_acq*delta_r)); %(2*pi)/gamma)*(1/(delta_r*t_acq))
BW_p = (6/0.001); %2n/tau
G_s = BW_p*(2*pi/(gamma*0.005));%slice thickness is 5mm 

%Generate the excitation(s)

%Generate the RF excitation waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulse(rfStepsE) = adipose_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla

%Generate the refocusing pule(s)
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
startR = (TE/2);

rfPulseR = adipose_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla
rfPulse(round((startR*10^5))+(rfStepsR)-1) = rfPulseR;

%Generate the readout + ADC
timestepsADCEvent = 576:623;
timestepsREADOUTEvent = 576:623;

adc(:,timestepsADCEvent) = 1:48;
gradAmp(1,timestepsREADOUTEvent) = G_r; %X gradients in Tesla per meter


% Generate the phase encoding gradient
timestepsPHASEPrephaseEvent = 526:575;
gradAmp(2,timestepsPHASEPrephaseEvent) = -(G_r/50)*24;  %Y gradients in Tesla per meter


%Slice refocusing gradient and read prephase
timestepsREADOUTPrephaseEvent = 526:575;

gradAmp(1,timestepsREADOUTPrephaseEvent) =  -(G_r/50)*24; %X gradients in Tesla per meter


%dephase and rephase grad
gradAmp(3,round(1:(pulsedurE/(dt)))) =  G_s; %Z gradients in Tesla per meter
gradAmp(3,round((pulsedurE/(dt))+1):(150)) =  -G_s; %Z gradients in Tesla per meter
gradAmp(3,round((startR*10^5))+(rfStepsR)-1) =  G_s; %Z gradients in Tesla per meter
gradAmp(3,250:299) =  -G_s; %Z gradients in Tesla per meter
gradAmp(3,400:449) =  -G_s; %Z gradients in Tesla per meter
% 

gradAmp(1,round(1:(pulsedurE/(dt)))) =  G_s; %Z gradients in Tesla per meter
gradAmp(1,round((pulsedurE/(dt))+1):(150)) =  -G_s; %Z gradients in Tesla per meter
gradAmp(1,round((startR*10^5))+(rfStepsR)-1) =  G_s; %Z gradients in Tesla per meter
gradAmp(1,250:299) =  -G_s; %Z gradients in Tesla per meter
gradAmp(1,400:449) =  -G_s; %Z gradients in Tesla per meter


%Plot the first TR
figure
subplot(5,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('TR Cycle: RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(5,1,2); plot(time,gradAmp(1,:),'r-','LineWidth',2);title('TR Cycle: Read Gradient'); ylim([-0.06 0.05]);
xlabel('time (s)'), ylabel('G_{r}(T/m)');grid on;
subplot(5,1,3); plot(time,gradAmp(2,:),'g-','LineWidth',2);title('TR Cycle: Phase Gradient');ylim([-0.06 0.05]);
xlabel('time (s)'), ylabel('G_{p}(T/m)');grid on;
subplot(5,1,4); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('TR Cycle: Slice Select Gradient');ylim([-0.06 0.05])
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

subplot(5,1,5); plot(time,adc,'k-','LineWidth',2);title('ADC_{x} Event');ylim([-20 60])
xlabel('time (s)'), ylabel('ADC Event Count');grid on;

%% Perform seuquence
kSpace  = zeros(size(T1,1),size(T1,2));
durPE = 0.5e-3;
gradPEmax = pi * (ySteps-1)/(gamma*FOV_r) * 1/durPE;

tic
for trIndex=1:size(T1,1)
    disp(trIndex);  
    
    %update the phase endcoding gradient
    gradAmp(2,timestepsREADOUTPrephaseEvent) =  (1-(trIndex-1)/((ySteps-1)/2)) * gradPEmax; %Y gradients in Tesla per meter
    k=1:xSteps;  
    j=1:ySteps;

        for i=1:zSteps            
            r = reshape(pos(:,k,j,i),[3,xSteps*ySteps]);
            mT = zeros(xSteps*ySteps,1);
            mZ = ones(xSteps*ySteps,1);
            
            dB0 = (gradAmp(:,1)'*r)'; 
            [mT,mZ] =  bloch(dt, dB0,rfPulse(1), reshape(T1(k,j),[xSteps*ySteps,1]),...
                                                 reshape(T2(k,j),[xSteps*ySteps,1]), mT, mZ);  
                
            for t=2:nTimeSteps 
                dB0 = (gradAmp(:,t)'*r)';
                [mT,mZ] =  bloch(dt, dB0,rfPulse(t), reshape(T1(k,j),[xSteps*ySteps,1]),...
                                                     reshape(T2(k,j),[xSteps*ySteps,1]), mT, mZ);
                
                if(adc(1,t)>0)
                    PDt = PD(k,j);
                    kSpace(round(adc(t)),trIndex) = kSpace(round(adc(t)),trIndex)+sum(mT.*reshape(PDt(k,j),[xSteps*ySteps,1]));
                end
                
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
