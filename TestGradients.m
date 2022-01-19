clear all
clc

LI = load('LineImage.mat');
SE = load('SpinEchoImage.mat');


figure;

subplot(2,2,1);imagesc(abs(SE.Ispace));

subplot(2,2,2);imagesc(abs(LI.Ispace));
title('Line Image');


subplot(2,2,4);imagesc(abs(SE.Ispace)-abs(LI.Ispace));
title('Test');