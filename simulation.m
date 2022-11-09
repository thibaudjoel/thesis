% vector of scaling parameters
alphas = (1:100)/100;
% edge probabilities for vertices of different communities in the unscaled connectivity matrix
q_s = [0.4 0.8];

% number of sampled adjacency matrices for each alpha
samplesize = 1000;

% size of community 1 and community 2 for n=50
[n1, n2] = deal(25, 25);
%compute mean error for n=50
mean_50_q_04 = sparsity_error(alphas, n1, n2, q_s(1), samplesize);
mean_50_q_08 = sparsity_error(alphas, n1, n2, q_s(2), samplesize);

% size of community 1 and community 2 for n=100
[n1, n2] = deal(50, 50);
%compute mean error for n=100
mean_100_q_04 = sparsity_error(alphas, n1, n2, q_s(1), samplesize);
mean_100_q_08 = sparsity_error(alphas, n1, n2, q_s(2), samplesize);

% size of community 1 and community 2 for n=150
[n1, n2] = deal(75, 75);
%compute mean error for n=150
mean_150_q_04 = sparsity_error(alphas, n1, n2, q_s(1), samplesize);
mean_150_q_08 = sparsity_error(alphas, n1, n2, q_s(2), samplesize);

% plot the mean errors for the different n's
plot(alphas, [mean_50_q_04, mean_100_q_04, mean_150_q_04,...
              mean_50_q_08, mean_100_q_08, mean_150_q_08],'LineWidth', 1)
 
legend('n=50, q=0.4', 'n=100, q=0.4', 'n=150, q=0.4',...
       'n=50, q=0.8', 'n=100, q=0.8', 'n=150, q=0.8');
xlabel('\alpha','fontsize',14);
ylabel('$\overline{L(\widehat{\Theta},\Theta)}$','fontsize',14,'interpreter','latex');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(gcf, 'mean_error.pdf');