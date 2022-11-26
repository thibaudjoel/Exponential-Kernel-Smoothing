rng(1000)
%% - Generate a function f from the library;
% number of coefficients included in random function
m=30; 
%generates random function f
f = randomfunction(m); 
%% - Generate data Y from the model for a standard Gaussian noise ε ∼ N (0, In)
% number of design points
n = 100;
% noise level
sigma = 1;
% hyperparameters/bandwidths
bandwidths = [0.001 0.01 0.05 0.1];
%% - Demonstrate (by figures) the behavior of the selected smoothing method fm with different
%   values of tuning parameter m ∈ M starting from undersmoothing (low complexity) and then
%   changing to oversmoothing (high complexity). Ensure that both effects meet.

%create design points
X = (0:n-1).'/n;
true_f = f(X);
% generate data from model
Y = f(X)+normrnd(0,sigma^2,[n,1]);
responses = zeros(n,length(bandwidths));
for i = 1:length(bandwidths)
    responses(:,i) = response(bandwidths(i),X,Y);
end
figure(1)
scatter(X,Y,30,[0.3 .5 .7],'filled','DisplayName','data')
hold on
plot(X,[f(X) responses],'LineWidth',1.5)
lgd = legend(['data' '$f$' strcat('$h$=',string(bandwidths))],'Interpreter','latex');
%title(lgd,'bandwidth h')
xlabel('$x$','Interpreter','latex','FontSize',16);
ylabel('$\tilde{f}_h(x)$','Interpreter','latex','FontSize',16);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'oversmoothing.pdf');

%% - Given f , compute by Monte-Carlo simulations and plot as function of m the mean fm =
%   E[fm] , the squared bias |f − fm|2, the variance tr Var(fm) , the squared risk Rm and the risk minimizer m∗

n_MC = 1000;
% hyperparameters/bandwidths
bandwidths = (1:50)/500;
X = (0:n-1).'/n;

Y_samples = zeros(n, n_MC);
for i = 1:n_MC
    Y_samples(:,i) = f(X)+normrnd(0,sigma^2,[n,1]);
end

mean_response = zeros(n, length(bandwidths));
squared_biases = zeros(1,length(bandwidths));
variances = zeros(1, length(bandwidths));
squared_risks = zeros(1, length(bandwidths));
unbiased_risk_estimates = zeros(n_MC, length(bandwidths));
lou_response_risk = zeros(n_MC,length(bandwidths));
unbiri_selectors = zeros(1,n_MC);
lou_selectors = zeros(1,n_MC);
counter = 0;
for h = bandwidths
    kernel_matrix=kern_mat(h,X);
    counter = counter+1;
    responses = zeros(n, n_MC);
    lou_responses = zeros(n, n_MC);
    for i = 1:n_MC
        responses(:,i) = response(h,X,Y_samples(:,i));
        lou_responses(:,i) = lou_response(h,X,Y_samples(:,i));
        unbiased_risk_estimates(i,counter) = norm(responses(:,i)-Y_samples(:,i))^2+2*sigma^2*trace(kernel_matrix) ;
        lou_response_risk(i,counter) = norm(lou_responses(:,i)-Y_samples(:,i))^2;
    end
    mean_response(:,counter) = mean(responses, 2);
    squared_biases(counter) = norm(f(X)-mean_response(:,counter))^2;
    variances(counter) = sigma^2*trace(kernel_matrix*kernel_matrix');
    squared_risks(counter) = squared_biases(counter)+variances(counter);
end
[minimal_risk, minimizing_bandwidth_index] = min(squared_risks);
risk_minimizer = bandwidths(minimizing_bandwidth_index);
[minimal_unbirisk, minimizing_unbibandwidth_index] = min(unbiased_risk_estimates,[],2);
[minimal_lourisk, minimizing_loubandwidth_index] = min(lou_response_risk,[],2);

squared_losses_unbiri = zeros(1,n_MC);
squared_losses_lou = zeros(1,n_MC);
for i = 1:n_MC
    unbiri_selectors(i) = bandwidths(minimizing_unbibandwidth_index(i));
    lou_selectors(i) = bandwidths(minimizing_loubandwidth_index(i));
    squared_losses_unbiri(i) = norm(response(unbiri_selectors(i),X,Y_samples(:,i))-f(X))^2;
    squared_losses_lou(i) = norm(response(lou_selectors(i),X,Y_samples(:,i))-f(X))^2;
end
figure(31)
plot(X, [f(X) mean_response(:,[10 20 30 40])],'LineWidth',1.5)
legend(['$f$' strcat('$h$=',string(bandwidths([10 20 30 40])))],'Interpreter','latex');
xlabel('$x$','Interpreter','latex','FontSize',16);
ylabel('$f_h(x)$','Interpreter','latex','FontSize',16);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'expresponse.pdf');
figure(40)
plot(bandwidths, [squared_biases;variances;squared_risks],'LineWidth',1.5)
legend('bias', 'variance', 'risk','Interpreter','latex')
xlabel('$h$','Interpreter','latex','FontSize',16);
xline(risk_minimizer,'--r','DisplayName','$h^*$')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'biasvariance.pdf');

     
%% - For each generated sample Y(k), implement the UnbiRiEs selector m^ and plot the obtained estimate fb against data,
%   true function, and oracle estimate fe m∗ ; see Figure 8.1;

%fix sample 
k=150;
Y_k = Y_samples(:,k);
unbiri_response = response(unbiri_selectors(k),X,Y_k);
lou_response = response(lou_selectors(k),X,Y_k);
oracle_response = response(risk_minimizer,X,Y_k);
figure(33);
scatter(X,Y_k,30,[0.3 .5 .7],'filled','DisplayName','data')
hold on
plot(X,f(X),'LineStyle','-','LineWidth',1.5);
plot(X,oracle_response,'LineStyle','--','LineWidth',1.5);
plot(X,unbiri_response,'LineStyle','-.','LineWidth',1.5);
plot(X,lou_response,'LineStyle',':','LineWidth',1.5);
legend('data','$f$','$\tilde{f}_{h^{*}}$','$\tilde{f}_{\hat{h}_{LO}}$','$\tilde{f}_{\hat{h}_{UR}}$','Interpreter','latex')
xlabel('$x$','Interpreter','latex','FontSize',16);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'datadrivestim.pdf');
%% - Plot the histogram of the selector m^ ; cf. Figure 8.4;
edges = linspace(0.002,0.1,50);
figure(7)
nbins=50;
histogram(unbiri_selectors,edges);
ylabel('count');
xlabel('$\hat{h}_{UR}$','Interpreter','latex','FontSize',16);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'histunbiri.pdf');

figure(8)
histogram(lou_selectors,edges);
ylabel('count');
xlabel('$\hat{h}_{LO}$','Interpreter','latex','FontSize',16);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'histlou.pdf');
%% - Illustrate (by plotting the obtained estimates) the variability of f_m∗ and of f^ using many Monte-Carlo realizations
% #MC realizations
l=20;
oracle_opt_responses = zeros(n,l);
lou_opt_responses = zeros(n,l);
unbiri_opt_responses = zeros(n,l);
for i=1:l
    oracle_opt_responses(:,i) = response(risk_minimizer,X,Y_samples(:,i));
    lou_opt_responses(:,i) = response(lou_selectors(i),X,Y_samples(:,i));
    unbiri_opt_responses(:,i) = response(unbiri_selectors(i),X,Y_samples(:,i));
end

%% - Plot the histogram of the squared loss |fb=^ − f|^2 against the oracle risk
figure(13)
edges_2 = linspace(0,100,30);
histogram(squared_losses_unbiri,edges_2);
xlabel('$\|\tilde{f}_{\hat{h}_{UR}}-f\|^2$','Interpreter','latex','FontSize',16);
ylabel('count');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'histlossunbiri.pdf');

figure(14)
histogram(squared_losses_lou,edges_2);
xlabel('$\|\tilde{f}_{\hat{h}_{LO}}-f\|^2$','Interpreter','latex','FontSize',16);
ylabel('count');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf,'histlosslou.pdf');

for i=1:l
    figure(i+100)
    scatter(X,Y_samples(:,i),30,[0.3 .5 .7],'filled','DisplayName','data')
    hold on
    plot(X,f(X),'LineStyle','-','LineWidth',1.5);
    plot(X,oracle_opt_responses(:,i),'LineStyle','--','LineWidth',1.5);
    plot(X,unbiri_opt_responses(:,i),'LineStyle','-.','LineWidth',1.5);
    plot(X,lou_opt_responses(:,i),'LineStyle',':','LineWidth',1.5);
    hold on

    xlabel('$x$','Interpreter','latex','FontSize',16);
    legend('data','$f$','$\tilde{f}_{h^*}$','$\tilde{f}_{\hat{h}_{UR}}$','$\tilde{f}_{\hat{h}_{LO}}$','Interpreter','latex')
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    saveas(gcf,strcat(string(i+100),'.pdf'));
end