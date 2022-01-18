%% Exercise 5.2
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


%% Question A

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

%Slice refocusing gradient and read prephase
timestepsREADOUTPrephaseEvent = 526:575;

gradAmp(1,timestepsREADOUTPrephaseEvent) =  -(G_r/50)*24; %X gradients in Tesla per meter
gradAmp(2,timestepsREADOUTPrephaseEvent) = -(G_r/50)*24;  %Y gradients in Tesla per meter
gradAmp(3,round(1:(pulsedurE/(dt)))) =  G_s; %Z gradients in Tesla per meter
gradAmp(3,round((pulsedurE/(dt))+1):(150)) =  -G_s; %Z gradients in Tesla per meter

gradAmp(3,round((startR*10^5))+(rfStepsR)-1) =  G_s; %Z gradients in Tesla per meter
gradAmp(3,250:299) =  -G_s; %Z gradients in Tesla per meter
gradAmp(3,400:449) =  -G_s; %Z gradients in Tesla per meter

%Plot the first TR
figure
subplot(5,1,1); plot(time,rfPulse,'k-','LineWidth',2);title('TR Cycle: RF Pulse'); 
xlabel('time (s)'), ylabel('|B_{1}^{+}| (T)');grid on;

subplot(5,1,2); plot(time,gradAmp(1,:),'r-','LineWidth',2);title('TR Cycle: Read Gradient'); ylim([-0.03 0.02]);
xlabel('time (s)'), ylabel('G_{r}(T/m)');grid on;
subplot(5,1,3); plot(time,gradAmp(2,:),'g-','LineWidth',2);title('TR Cycle: Phase Gradient');ylim([-0.03 0.02]);
xlabel('time (s)'), ylabel('G_{p}(T/m)');grid on;
subplot(5,1,4); plot(time,gradAmp(3,:),'b-','LineWidth',2);title('TR Cycle: Slice Select Gradient');ylim([-0.06 0.1])
xlabel('time (s)'), ylabel('G_{z}(T/m)');grid on;

subplot(5,1,5); plot(time,adc,'k-','LineWidth',2);title('ADC_{y} Event');ylim([-20 60])
xlabel('time (s)'), ylabel('ADC Event Count');grid on;

%% Question B:

kSpace  = zeros(size(T1,1),size(T1,2));

tic
for trIndex=1:size(T1,2)

    disp(trIndex);    

    %Update phase encoding gradients
    gradAmp(2,timestepsREADOUTPrephaseEvent) =  (-(G_r/50)*24)+(((G_r/50)*24)/48)*(trIndex-1); %Y gradients in Tesla per meter


    k=1:xSteps;
    for j=1:ySteps
        for i=1:zSteps
            r = pos(:,k,j,i);
            mT = zeros(xSteps,1);
            mZ = ones(xSteps,1);

            %dB0 = ...
            dB0 = gradAmp(:,1)'*r; 

%             dB0 = dot(gradAmp(:,1)',r); 
            [mT,mZ] =  bloch(dt, dB0,rfPulse(1), T1(k,j),T2(k,j), mT, mZ);

            for t=2:nTimeSteps %i starts at 2

%                 dB0 = dot(gradAmp(:,t),r(:,j)); 
                dB0 = gradAmp(:,t)'*r; 
                [mT,mZ] =  bloch(dt, dB0,rfPulse(t), T1(k,j),T2(k,j), mT, mZ);  
                
                
                if(adc(t)>0)
                    kSpace(round(adc(t)),trIndex) = kSpace(round(adc(t)),trIndex) + (sum(mT)*PD(k,j));
%                     fprintf('index %d,%d\n',round(adc(t)),trIndex);
                end

            end            


        end
    end

end
toc

%% Question C: Reconstruct the image

figure;
Ispace = fftshift((ifft2(kSpace)));
subplot(2,1,1);imagesc(abs(Ispace));
title('Magnitude Image');


subplot(2,1,2);imagesc(log(abs((kSpace))),'CDataMapping','scaled'); 
title('Absolute kSpace');
saveas(gcf,'5p2')

