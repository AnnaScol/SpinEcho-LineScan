clear all
clc
%%
LI_1D = load('LineImage1D.mat');
LI = load('LineImage.mat');
% SE = load('SpinEchoImage.mat');


figure;

line2D = LI.Ispace(:,25)';


subplot(2,2,1);plot(normalize(abs(line2D)));
title('2D Line Image');
subplot(2,2,2);plot(normalize(abs(LI_1D.Ispace))');
title('1D Line Image');


subplot(2,2,4);plot(normalize(abs(line2D))-normalize(abs(LI_1D.Ispace))');
title('Difference');


figure;

subplot(2,1,1);imagesc(normalize(abs(line2D)));
title('2D Line Image');
subplot(2,1,2);imagesc(normalize(abs(LI_1D.Ispace)'));
title('1D Line Image');

%% kSpace comparison
k1D = fftshift(fft(fftshift(LI_1D.Ispace)));
kLine = fftshift(fft(fftshift(line2D)));

figure;
subplot(2,1,1);imagesc(log(abs(k1D)));
title('1D Line Image');

subplot(2,1,2);imagesc(log(abs(kLine)));
title('From 2D Line Image');
