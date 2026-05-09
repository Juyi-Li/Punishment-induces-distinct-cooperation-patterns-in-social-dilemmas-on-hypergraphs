%clc; clear; close all;

% ---------------------------
% Data
% ---------------------------
alphas = 0:0.01:1;
betas  = 0:0.01:1;

Z = CE_domain';   


% Check whether the dimensions are consistent
if size(Z,1) ~= length(alphas) || size(Z,2) ~= length(betas)
    error('The size of CE_domain should be length(alphas) x length(betas).');
end

% ---------------------------
% Take the value corresponding to (alpha,beta) = (0,0)
% Since alpha = 0 corresponds to the first row and beta = 0 corresponds to the first column
% ---------------------------
Z00 = 0.001;

% Global color range
cmin = min(Z(:));
cmax = max(Z(:));

if cmax == cmin
    error('All values in Z are identical, so an effective colorbar cannot be generated.');
end

% ---------------------------
% Customized colormap
% Make Z00 correspond to a nearly white color
% Color scheme: sand-like color -> near white -> dark teal
% ---------------------------

lowColor = [0.96 0.95 0.93];
%lowColor   = [0.86 0.74 0.46];   % Sand yellow
whiteColor = [0.96 0.95 0.93];   % Near white
highColor  = [0.00 0.30 0.27];   % Dark teal

% Position of Z00 in the colorbar, ranging from 0 to 1
t0 = (Z00 - cmin) / (cmax - cmin);
t0 = max(0, min(1, t0));

n  = 256;
n1 = max(2, round((n-1)*t0));   % Number of colors below white
n2 = max(2, n - n1);            % Number of colors above white

% Lower segment: sand-like color -> white
cmap1 = [linspace(lowColor(1), whiteColor(1), n1)', ...
         linspace(lowColor(2), whiteColor(2), n1)', ...
         linspace(lowColor(3), whiteColor(3), n1)'];

% Upper segment: white -> dark teal
cmap2 = [linspace(whiteColor(1), highColor(1), n2)', ...
         linspace(whiteColor(2), highColor(2), n2)', ...
         linspace(whiteColor(3), highColor(3), n2)'];

cmap = [cmap1; cmap2];

% ---------------------------
% Plot
% ---------------------------
figure('Color','w');

% Note: the horizontal axis is beta, and the vertical axis is alpha
imagesc(betas, alphas, Z);
axis xy;
hold on;

colormap(cmap);
caxis([cmin cmax]);

% Colorbar
cb = colorbar;
cb.LineWidth = 1.2;
cb.FontSize = 14;
cb.Label.String = 'CE\_domain';
cb.Label.FontSize = 16;

% ---------------------------
% Mark the point (alpha,beta) = (0,0)
% The horizontal coordinate is beta, and the vertical coordinate is alpha
% ---------------------------
% plot(0, 0, 'ko', ...
%     'MarkerSize', 7, ...
%     'LineWidth', 1.5, ...
%     'MarkerFaceColor', 'w');
% 
% text(0.03, 0.04, '(0,0)', ...
%     'FontSize', 12, ...
%     'Color', 'k', ...
%     'FontWeight', 'bold');
% 
% plot(betas, betas, 'k--', 'LineWidth', 2.2);

% ---------------------------
% Plot the contour line for Z = Z00
% In contour(X,Y,Z):
% X corresponds to the column variable, i.e., beta
% Y corresponds to the row variable, i.e., alpha
% ---------------------------

% contour(betas, alphas, Z, [Z00 Z00], ...
%     'k--', 'LineWidth', 2.2);

%mask = Z > 0;
%contour(betas, alphas, double(mask), [0.5 0.5], 'k--', 'LineWidth', 2.2);

% ---------------------------
% Axis settings
% ---------------------------
xlabel('\delta', 'FontSize', 20);
ylabel('\alpha', 'FontSize', 20);

set(gca, ...
    'FontSize', 16, ...
    'LineWidth', 1.8, ...
    'Box', 'on', ...
    'Layer', 'top');

pbaspect([1 1 1]);

xlim([min(betas), max(betas)]);
ylim([min(alphas), max(alphas)]);