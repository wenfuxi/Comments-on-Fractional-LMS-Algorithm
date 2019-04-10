% LMS and fractional LMS for System Identification

clear all
close all
clc

tic

num = 30000; % number of measurements

dim = 16; % dimension of the parameters

gain = 10;
w = gain + rand(dim,1)*1*gain;
% w = gain*randn(dim,1);

% w = [-2 3].';
% w = gain + rand(dim,1)*gain;
% w = ones(dim,1);
 
len = num + dim - 1; % number of input signals

input = randn(len,1); % input signal

x = zeros(dim,num); % signal in matrix form

for it = 1:num
    x(:,it) = input(it:it+dim-1,1);
end

noise = randn(1,num);

snr = 20; % in dB scale

d = w.'*x + 10^(-snr/20)*noise;

mu1 = 0.0005;
mu2 = 0.0005;
fp = 0.5; % fractional power

mu = 0.0005;

Rdx = mean(ones(dim,1)*d.*x,2);
Rxx = x*x'/num;
wot = Rxx\Rdx;
numt = 100000;

sigt = gain*(2*rand(1,numt)-1);
% sigt = gain + gain*ones(1,numt);

t1 = abs(mean(sigt.^0.5));
t2 = abs(mean(sigt.^1.5));

wo = mu*10^(-snr/10)*trace(Rxx)/2;

wo2 = mu1*(wo*(1+t1/gamma(1.5)) - t2/gamma(1.5));
% Eopt(1,it) = norm(wo - w)^2;
 
% wo2 = (wo*(1+t1/gamma(1.5))- t2/gamma(1.5));

Wn = gain + ones(dim,1)*1*gain;
e = zeros(1,num);
Err = zeros(1,num);

Wn2 = gain + ones(dim,1)*1*gain;
e2 = zeros(1,num);
Err2 = zeros(1,num);

WnA = gain + ones(dim,1)*1*gain;
eA = zeros(1,num);
ErrA = zeros(1,num);

for it = 1:num
    % LMS algorithm    
    e(1,it) = d(1,it) - Wn(:,it)'*x(:,it);
    Wn(:,it+1) = Wn(:,it) + mu*e(1,it)*x(:,it);
    Err(1,it) = norm(Wn(:,it+1) - w)^2;
    
    % Fractional LMS algorithm - REAL
%     Wn2(:,it) = real(Wn2(:,it));
    e2(1,it) = d(1,it) - Wn2(:,it)'*x(:,it);
    temF = Wn2(:,it).^(1-fp);
    Wn2(:,it+1) = Wn2(:,it) + mu1*e2(1,it)*x(:,it) + mu2*e2(1,it)*(x(:,it).*temF)/gamma(2-fp); 
    Err2(1,it) = norm(real(Wn2(:,it+1)) - w(:,1))^2;
    
        % Fractional LMS algorithm - ABS
    WnA(:,it) = abs(WnA(:,it));
    eA(1,it) = d(1,it) - WnA(:,it)'*x(:,it);
    temFA = WnA(:,it).^(1-fp);
    WnA(:,it+1) = WnA(:,it) + mu1*eA(1,it)*x(:,it) + mu2*eA(1,it)*(x(:,it).*temFA)/gamma(2-fp); 
    ErrA(1,it) = norm(real(WnA(:,it+1)) - w(:,1))^2;
    
    %Eopt(1,it) = norm(wo - w)^2;
    
    Eopt(1,it) = wo;
    Eopt2(1,it) = norm(wo2 - w)^2;
%     Eopt2(1,it) = wo2;
end

figure

px = 1:num;
semilogy(px,Err2,'-',px,ErrA,'-',px,Eopt2/2,'k-.',px,Err,px,Eopt,'k--','linewidth',2) 
legend('Fractional LMS: Real','Fractional LMS: Abs','Fractional LMS (theoretical)','LMS','LMS (theoretical)')

xlabel('iteration','fontsize',12)
ylabel('Mean squared deviation (MSD)','fontsize',12)
% title(['w = [',num2str(w.'),' ], \mu = \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(fp)],'fontsize',12)

figure

% subplot(121)
plot(real(Wn2.'),'linewidth',2,'DisplayName','Wn2'),grid
xlabel('iteration','fontsize',12)
ylabel('w(n)','fontsize',12)
% title(['w = [',num2str(w.'),' ]'],'fontsize',12)
axis([1,num,min(w)-1,max(w)+1])

% subplot(122)
% plot(imag(Wn2.'),'linewidth',2,'DisplayName','Wn2'),grid
% xlabel('iteration','fontsize',12)
% ylabel('Imaginary part of w(n)','fontsize',12)
% axis([1,num,min(w)-1,max(w)+1])
% title(['w = [',num2str(w.'),' ]'],'fontsize',12)

figure
 
% subplot(121)
plot(real(WnA.'),'linewidth',2,'DisplayName','Wn2'),grid
xlabel('iteration','fontsize',12)
ylabel('w(n)','fontsize',12)
%title(['w = [',num2str(w.'),' ],fractional LMS: Abs'],'fontsize',12)
axis([1,num,min(w)-1,max(w)+1])

% subplot(122)
% plot(imag(WnA.'),'linewidth',2,'DisplayName','Wn2'),grid
% xlabel('iteration','fontsize',12)
% ylabel('Imaginary part of w(n)','fontsize',12)
% axis([1,num,min(w)-1,max(w)+1])
% title(['w = [',num2str(w.'),' ]'],'fontsize',12)

% figure
% plot(real(temp.'),'linewidth',2),grid
% xlabel('iteration','fontsize',12)
% title(['w(n) = [',num2str(w.'),' ]'],'fontsize',12)
toc