% LMS and fractional LMS for System Identification

clear all
close all
clc

tic

num = 15000; % number of measurements

dim = 8; % dimension of the parameters

w = 10 + rand(dim,1)*10;
% w = randn(dim,1);

len = num + dim - 1; % number of input signals

input = randn(len,1); % input signal

x = zeros(dim,num); % signal in matrix form

for it = 1:num
    x(:,it) = input(it:it+dim-1,1);
end

noise = randn(1,num);

snr = 20; % in dB scale

d = w.'*x + 10^(-snr/10)*noise;

Rdx = mean(ones(dim,1)*d.*x,2);
Rxx = mean(x.*x,2);
wo = Rdx./Rxx;
% Eopt(1,it) = norm(wo - w)^2;

mu = 0.0025;

mu1 = 0.0005;
mu2 = mu1;

fpL = 0.2;
fpS = 0.2;
fpH = 1;
% fp = 0.5; % fractional power

Wn = randn(dim,1);
e = zeros(1,num);
Err = zeros(1,num);

Wn2 = randn(length(fpL:fpS:fpH),dim,2);
e2 = zeros(1,num);
Err2 = zeros(1,num);

for it = 1:num
    % LMS algorithm    
    e(1,it) = d(1,it) - Wn(:,it)'*x(:,it);
    Wn(:,it+1) = Wn(:,it) + mu*e(1,it)*x(:,it);
    Err(1,it) = norm(Wn(:,it+1) - w)^2;
    
    itp = 0;
    for fp = fpL:fpS:fpH
        itp = itp + 1;
        % Fractional LMS algorithm    
        e2(itp,it) = d(1,it) - Wn2(itp,:,it)*x(:,it);
        Wn2(itp,:,it) = abs(Wn2(itp,:,it));
        temF = Wn2(itp,:,it).^(1-fp);
        Wn2(itp,:,it+1) = Wn2(itp,:,it).' + mu1*e2(itp,it)*x(:,it) + mu2*e2(itp,it)* (x(:,it).*temF.')/gamma(2-fp);
        Err2(itp,it) = norm(Wn2(itp,:,it+1).' - w)^2;
    end
    
%     Eopt(1,it) = norm(wo - w)^2;
    end

% figure
px = 1:num;
% semilogy(px,Err,px,Err2,'-.','linewidth',1.5),grid
% legend(['LMS: \mu = ',num2str(mu)],['Fractional LMS: \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(fp)])
% 
% xlabel('iteration','fontsize',12)
% ylabel('Mean squared deviation (MSD)','fontsize',12)

figure
px = 1:num;
semilogy(px,Err2(1,:),'--',px,Err2(2,:),'-.',px,Err2(3,:),':',px,Err2(4,:),'-',px,Err2(5,:),'--','linewidth',1.5),
str1 = ['Fractional LMS: \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(0.2)];
str2 = ['Fractional LMS: \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(0.4)];
str3 = ['Fractional LMS: \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(0.6)];
str4 = ['Fractional LMS: \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(0.8)];
str5 = ['Fractional LMS: \mu_1 = \mu_2 = ',num2str(mu1),', v = ',num2str(1)];
legend(str1,str2,str3,str4,str5),grid
xlabel('iteration','fontsize',12)
ylabel('Mean squared deviation (MSD)','fontsize',12)

toc