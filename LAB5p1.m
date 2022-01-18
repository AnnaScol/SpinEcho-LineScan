%% Exercise 5.1
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-5; 
gamma = 2*pi*42.577*10^6;

te = [6 8 12 16 24 32 48 64 96 128 192 256]*10^-3; % te in ms

%Question A

%Allocate the memory needed
nTimeSteps  = 70000;
rfPulseE    = zeros(1,100); %excitation
rfPulseR    = zeros(1,100); %refcosuing
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,nTimeSteps); %variable to hold 10000 time points

tiSteps    = size(te,2);
nSpins     = 1000;
spins      = zeros(nSpins,1); %variable to hold nPosSteps positions allong the z direction
mFinalVect = zeros(tiSteps,2); %variable to hold the final magnetization calculated for each position

%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

%Generate off-resonances 
for i=1:nSpins
    spins(i) = 2.5*(rand-0.5)*10^-7;     
end

% relaxation times
T1 = 0.850; %s
T2 = 0.06; %s

%% Question A:

%Generate the RF excitation waveform
pulsedurE = 0.001; % duration of the RF in s
rfStepsE = round(1:(pulsedurE/(dt)));
rfPulseE(rfStepsE) = adipose_sinc_rf(length(rfStepsE),3,pi/2,dt); %B1+ in Tesla


%Generate the refocusing waveform
pulsedurR = 0.001; % duration of the RF in s
rfStepsR = round(1:(pulsedurR/(dt)));
rfPulseR(rfStepsR) = adipose_sinc_rf(length(rfStepsR),3,pi,dt); %B1+ in Tesla


%% Question B:


%Question B:

for measIndex = 1:size(te,2)

%Generate the sequence
    startE = 1;
    startR = (te(measIndex)/2) + 0.0005;
    endTime = te(measIndex)+(pulsedurE/2);
    
    rfPulse     = zeros(1,nTimeSteps); 
    rfPulse((startE-1)+(rfStepsE)) = rfPulseE;
    rfPulse(round((startR*10^5))+(rfStepsR)-1) = rfPulseR;
 
    mT = zeros(nSpins,1);
    mZ = ones(nSpins,1);
        
    %simulate the sequence
    spinIndex=1:nSpins;
    gradAtPosJ = spins(spinIndex);
    [mT,mZ] =  bloch(dt, gradAtPosJ,rfPulse(1), T1, T2, mT, mZ);
    for i=2:(endTime*10^5) %i starts at 2
        [mT,mZ] =  bloch(dt, gradAtPosJ,rfPulse(i), T1, T2, mT, mZ);
    end
    
    
    mFinalVect(measIndex,1) = mean(mT);
    mFinalVect(measIndex,2) = mean(mZ);
end

%% Question C:
figure;
subplot(2,1,1);plot(te,abs(mFinalVect(:,1)),'o:','LineWidth',2);
xlabel('TE (s)'), ylabel('|M_{xy}|');
title('T_2-relaxation');
subplot(2,1,2);plot(te,angle(mFinalVect(:,1)),'o:','LineWidth',2);
xlabel('TE (s)'), ylabel('\angleM_{xy} ')
title('T_2-relaxation')

fo = fitoptions('Method','NonlinearLeastSquares');
ft = fittype('a*exp(-x/b)','options',fo);
[curve,gof,output] = fit(te',abs(mFinalVect(:,1)),ft);

figure;plot(te',abs(mFinalVect(:,1)),'o','LineWidth',2);
hold on;plot(te',abs(curve.a*exp(-te'/curve.b)),'-','LineWidth',2);
xlabel('TE (s)'), ylabel('|M_{xy}|')
title(['T_2-fit = ' num2str(curve.b*1000) 'ms'])
legend('S','T_2^*-fit','Location','best');legend boxoff;box on;


