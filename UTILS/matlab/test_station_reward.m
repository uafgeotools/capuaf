clear all
close all

Nmax = 50;
N = [1:Nmax];
max_penalty = 1.5;
min_penalty = 0.5;
N_threshold = 8;
font_size = 13;

cmap = jet(10);
ichosen = 7; % plot the chosen one
cmap(ichosen,:) = cmap(ichosen,:)*0;
mark_vec = {'--' '--' '--' '--' '--' '--' '-' '--' '--' '--'};
figure; hold on

for ii = 1:1:10
    %k = (1 - (2/pi)*atan(N/ii));
    %k = (1 - (2/pi)*atan(N/ii))*ii;
    %k = (1 - (2/pi)*atan(N))*ii;
    k(ii,:) = (exp(-N/ii)* (1+min_penalty))+ min_penalty;
    plot(N,k(ii,:),mark_vec{ii},'LineWidth',3,'Color', cmap(ii, :));
end

set(gca,'fontsize',font_size) 
xlabel('Number of stations (Ns)');
ylabel('Station reward factor (h)');
title('h(Ns) = 0.5 + 1.5exp(-Ns/C)');
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
set(v,'string','C','fontsize',font_size);

% save figure as eps
pdir = '/home/vipul/manuscripts/2016/nehrp/figures/cap/misfit_tests/';
ftag = 'Nstn_reward';
orient portrait; print(gcf,'-depsc',sprintf('%s%s',pdir,ftag));