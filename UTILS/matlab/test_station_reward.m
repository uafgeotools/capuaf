clear all
close all

Nmax = 50;
N = [1:Nmax];
max_penalty = 1.5;
min_penalty = 0.5;
N_threshold = 8;
font_size = 13;

figure; hold on

for ii = 1:1:10
    %k = (1 - (2/pi)*atan(N/ii));
    %k = (1 - (2/pi)*atan(N/ii))*ii;
    %k = (1 - (2/pi)*atan(N))*ii;
    k(ii,:) = (exp(-N/ii)* (1+min_penalty))+ min_penalty;
end
plot(N,k,'--','LineWidth',3);
set(gca,'fontsize',font_size) 
xlabel('Number of stations (Nstn)');
ylabel('Station penalty factor (k)');
title('k = (exp(-Nstn/CONST)*1.5)+0.5');
x1 = N_threshold*ones(1,100);
y1 = linspace(min_penalty,max_penalty,100);
y2 = ones(1,100);
x2 =linspace(0,Nmax,100);
hold on
plot(x1,y1);
plot(x2,y2);
axis ([0 30 min_penalty max_penalty]);

h = legend('1','2','3','4','5','6','7','8','9','10');
v = get(h,'title');
set(v,'string','CONST','fontsize',font_size);