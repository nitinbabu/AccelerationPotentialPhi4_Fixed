%% Post Pro


normxforce =    1000*9.81*0.25*Amp;
normzforce =    1000*9.81*0.5*Amp;
normymoment =  1000*9.81*0.5*0.25*Amp;

figure(65);
sgtitle('Forces on Surface Piercing body B=0.5m d = 0.25m, \lambda = 2m')
subplot(3,1,1)
title('Horizontal Force')
plot(timeaxis(1:end), Fx,'r-^');
hold on
legend('Fx')
grid minor
ylabel('Force-x')
xlabel('t(s)')
subplot(3,1,2)

title('Vertical Force')
plot(timeaxis(1:end), Fz,'b-^');
legend('Fz')
hold on
grid minor
ylabel('Force-z')
xlabel('t(s)')
%

subplot(3,1,3)

title('Moment')
plot(timeaxis(1:end-1), My(1:end-1),'b-^');
legend('My')
hold on
grid minor
ylabel('Moment-y')
xlabel('t(s)')



Fs = 1/dt;
figure(66);
subplot(1,3,2)
% Step 2: Calculate the FFT
fft_dataz = fft(Fz(400:end-1)./normzforce);
% Step 3: Compute the amplitudes and frequencies
N = length(fft_dataz);
amplitudes = 2/N * abs(fft_dataz(1:N/2+1)); % Compute the one-sided amplitudes
frequencies = (0:N/2) * Fs / N; % Compute the corresponding frequencies

% Plot the amplitude spectrum
% title('')
plot(frequencies(2:end), amplitudes(2:end));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum');

legend('Fz')
% Note: The amplitudes are scaled by 2/N to account for the one-sided amplitude spectrum.
% The amplitudes are doubled to represent the energy in both the positive and negative frequency components.

% You can also find the dominant frequency (peak frequency) using the following:
[max_amplitude, index] = max(amplitudes);
dominant_frequency = frequencies(index);
disp(['Dominant frequency: ', num2str(dominant_frequency), ' Hz']);

subplot(1,3,1)  %% X force

% Step 2: Calculate the FFT
fft_dataX = fft(Fx(400:end-1)./normxforce);
% Step 3: Compute the amplitudes and frequencies
N = length(fft_dataX);
amplitudes = 2/N * abs(fft_dataX(1:N/2+1)); % Compute the one-sided amplitudes
frequencies = (0:N/2) * Fs / N; % Compute the corresponding frequencies

% Plot the amplitude spectrum

plot(frequencies(2:end), amplitudes(2:end));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum');
legend('Fx')
% Note: The amplitudes are scaled by 2/N to account for the one-sided amplitude spectrum.
% The amplitudes are doubled to represent the energy in both the positive and negative frequency components.

% You can also find the dominant frequency (peak frequency) using the following:
[max_amplitude, index] = max(amplitudes);
dominant_frequency = frequencies(index);
disp(['Dominant frequency: ', num2str(dominant_frequency), ' Hz']);


subplot(1,3,3)
% Step 2: Calculate the FFT
fft_dataX = fft(My(400:end-1)./normymoment);
% Step 3: Compute the amplitudes and frequencies
N = length(fft_dataX);
amplitudes = 2/N * abs(fft_dataX(1:N/2+1)); % Compute the one-sided amplitudes
frequencies = (0:N/2) * Fs / N; % Compute the corresponding frequencies

% Plot the amplitude spectrum

plot(frequencies(2:end), amplitudes(2:end));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum');
legend('Fx')
% Note: The amplitudes are scaled by 2/N to account for the one-sided amplitude spectrum.
% The amplitudes are doubled to represent the energy in both the positive and negative frequency components.

% You can also find the dominant frequency (peak frequency) using the following:
[max_amplitude, index] = max(amplitudes);
dominant_frequency = frequencies(index);
disp(['Dominant frequency: ', num2str(dominant_frequency), ' Hz']);

%%


%%
figure;
quiver(bodynodexbody, bodynodeybody, NXbdiff, NYbdiff)
hold on
plot(bodynodexbody, bodynodeybody,'o-r')
axis equal
grid minor