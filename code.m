% %*******************************************************************************
% %********************************* TIME DOMAIN *******************************
% %*******************************************************************************
% 
clc ;
clear all;
% 
% % ================================
% % ======== reading audios ========
% % ================================
% 
[x1,fs1]=audioread('Team2_ speechsignal_1.wav');
[x2,fs2]=audioread('Team2_ speechsignal_2.wav');
[x3,fs3]=audioread('Team2_ speechsignal_3.wav');

L1 = length(x1);  
L2 = length(x2);  
L3 = length(x3);  
x1=([x1' zeros(1,L3-L1)])';
x2=([x2' zeros(1,L3-L2)])';
x3=([x3' zeros(1,L3-L3)])';

x1=interp(x1,20);
x2=interp(x2,20);
x3=interp(x3,20);

fs1=fs1*20;
TS1 = 1/fs1;                            
N1 = [0:(length(x1)-1)];
L1 = length(x1);                        
f1 = [-L1/2:L1/2-1]*(fs1/L1);           

fs2=fs2*20;
TS2 = 1/fs2;                            
N2 = [0:(length(x2)-1)];
L2 = length(x2);                        
f2 = [-L2/2:L2/2-1]*(fs2/L2);           

fs3=fs3*20;
TS3 = 1/fs3;                            
N3 = [0:(length(x3)-1)];
L3 = length(x3);                        
f3 = [-L3/2:L3/2-1]*(fs3/L3);           

% % ================================
% % ========== Modulating ==========
% % ================================
Fc1=50000;
carrier1 = cos(2*pi*Fc1*TS1*N1);
x1_modulated = x1'.* carrier1;
Fc2=200000;
carrier2 = cos(2*pi*Fc2*TS2*N2);
x2_modulated = x2'.*carrier2;
Fc3=200000;
carrier3 = sin(2*pi*Fc3*TS3*N3);
x3_modulated = x3'.*carrier3;

% figure();
% subplot(3,1,1);
% plot(f1,abs(fftshift(fft(x1_modulated))));
% title('X1 Modulated')
% subplot(3,1,2);
% plot(f2,abs(fftshift(fft(x2_modulated))));
% title('X2 Modulated')
% subplot(3,1,3);
% plot(f3,abs(fftshift(fft(x3_modulated))));
% title('X3 Modulated')


 
% %*******************************************************************************
% %******************************* COMBINED SIGNAL *******************************
% %******************************************************************************* 
signal = x1_modulated + x2_modulated + x3_modulated ;
figure();
plot(signal);
title('Time Domain')
figure();
plot(f3,abs(fftshift(fft(signal))));
title('Frequency Domain')

[C1,R1]=size(x1);
[C2,R2]=size(x2);
[C3,R3]=size(x3);

% %*******************************************************************************
% %******************************* Demodulate ************************************
% %*******************************************************************************

for i = [0 pi/18 pi/6 pi/2]
    carrier1 = cos(2*pi*Fc1*TS1*N1 + i);
    x1_new = signal(1:R1, 1:C1) .* carrier1;
    lp = designfilt('lowpassfir', 'FilterOrder',128, 'CutoffFrequency',22*10^3, 'SampleRate', fs1);%lowpassfilter
    x1_new=filter(lp,x1_new);
    OUTPUT1=downsample(x1_new,20);
    audiowrite(strcat('audio1Out',int2str(i*180/pi),'.wav'),OUTPUT1,fs1/20);


    carrier2 = cos(2*pi*Fc2*TS2*N2 + i);
    x2_new = signal(1:R2, 1:C2) .* carrier2;
    lp = designfilt('lowpassfir', 'FilterOrder',64, 'CutoffFrequency',44*10^3, 'SampleRate', fs2);%lowpassfilter
    x2_new=filter(lp,x2_new);
    OUTPUT2=downsample(x2_new,20);
    audiowrite(strcat('audio2Out',int2str(i*180/pi),'.wav'),OUTPUT2,fs2/20);

    carrier3 = sin(2*pi*Fc3*TS3*N3 + i);
    x3_new = signal(1:R3, 1:C3) .* carrier3;
    lp = designfilt('lowpassfir', 'FilterOrder',64, 'CutoffFrequency',44*10^3, 'SampleRate', fs3);%lowpassfilter
    x3_new=filter(lp,x3_new);
    OUTPUT3=downsample(x3_new,20);
    audiowrite(strcat('audio3Out',int2str(i*180/pi),'.wav'),OUTPUT3,fs3/20);
    
%     figure();
%     subplot(3,1,1);
%     plot(f1,abs(fftshift(fft(x1_new))));
%     title(strcat('X1 De-Modulated ',int2str(i*180/pi)))
%     subplot(3,1,2);
%     plot(f2,abs(fftshift(fft(x2_new))));
%     title(strcat('X2 De-Modulated ',int2str(i*180/pi)))
%     subplot(3,1,3);
%     plot(f3,abs(fftshift(fft(x3_new))));
%     title(strcat('X3 De-Modulated ',int2str(i*180/pi)))
    
end;
