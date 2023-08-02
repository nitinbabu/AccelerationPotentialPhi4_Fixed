function [Seta,Amp,f,k,eph]=spec_JONSWAP2(Hs,Ts,h)

g = 9.81; 
N = 101;
gamma = 3.3;        %peakedness param

m0 = (Hs/4)^2;
HMEAN = 2.5*sqrt(m0);

Tp = Ts/(1-0.132*(gamma+0.2)^-0.559); % equation 6 
sigma1 = 0.2189; % enhancement factor
fp = 1/Tp; fs = 1/Ts;  ws = 2*pi*fs;

sigmaA = 0.07; sigmaB = 0.09;

fmin = 0.5;  % minimum cutoff frequency
fmax = 2; % maximum cutoff frequency

f = zeros(1,N); f(1) = fmin; f(N) = fmax; Seta = zeros(1,N); Amp = zeros(1,N); w = zeros(1,N); dfk = zeros(1,N);
lambda =  zeros(1,N);
for n = 1:N-1
    dfk(n) = ( (fmax/fmin)^(1/ (N-1) )  - 1)*f(n);
    f(n+1) = f(n) + dfk(n);
end

for n = 1:N
    if f(n) <= fp 
        sigma0 = sigmaA;
    else
        sigma0 = sigmaB;
    end
    Seta(n) = sigma1*(Hs/Tp/Tp)^2*f(n)^-5*exp( -1.25*(Tp*f(n))^-4 )*gamma^exp(-(Tp*f(n)-1)^2/(2*sigma0^2) );
    Amp(n) = sqrt(2*Seta(n)*dfk(n));
    w(n) = 2*pi*f(n);
    [lam]=dispersioncalc(w(n),h);
    lambda(n) = lam;
    k(n) = 2*pi/lambda(n);
    
end
wave_ph=2*pi*rand(N,1)'; eph = wave_ph;


figure(10)
sgtitle('Hs = 0.03m  Ts = 1s')
subplot(2,2,1);
plot(f,Seta.*1e4)
xlabel('f (Hz)')
ylabel('S_\eta(f) cm^2/Hz')
grid on
title('JONSWAP Spectrum \gamma = 3.3')
grid minor

subplot(2,2,2);
plot(f,Amp)
xlabel('f (Hz)')
ylabel('Amp (m)')
grid on
title('Wave Height Spectrum')
grid minor
hold on
% plot(fs,Hs/2,'d')

subplot(2,2,3);
plot(f,lambda)
xlabel('f (Hz)')
ylabel('\lambda (m)')
grid on
hold on
[Ls]=dispersioncalc(ws,h);
plot(fs,Ls,'d')
legend('\lambda','\lambda_S')
title('\lambda vs f')
grid minor

subplot(2,2,4);
plot(f,h./lambda)
xlabel('f (Hz)')
ylabel('h / \lambda ')
grid on
hold on
plot(f,0.2+zeros(1,N),'-k')
legend('h/L','h/L = 0.2')
title('d/L ratio')
grid minor

%% Time Series
figure(5)
t = 0:0.02:600;
eta = 0;
L = 6*Ls;
for i=1:N
   eta = eta + Amp(i)*cos(k(i)*L/2-w(i).*t + eph(i));
end
plot(t,eta)
grid minor
xlabel('time (s)')
ylabel('\eta (m)')

%% Elevation
figure(6)
hold on

xf = 0:L/120:L;
tend = 50;
etax = 0;
for i=1:N-1
   etax = etax + Amp(i)*cos(k(i).*xf-w(i)*tend + eph(i));
end
plot(xf,etax)
grid minor

xlabel('xf ')
ylabel('\eta (m)')


