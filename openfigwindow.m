function [fig] = openfigwindow(UIAxes,visible)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

set(groot, 'DefaultTextInterpreter', 'LaTeX', ...
           'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
           'DefaultAxesFontName', 'LaTeX', ...
           'DefaultLegendInterpreter', 'LaTeX', ...
           'defaultFigureColor','w');

fig = figure;
if nargin == 2
    if visible == 'off'
        fig.Visible = 'off';
    end
end
figAxes = axes(fig);
set(figAxes,'FontSize',UIAxes.XLabel.FontSize);
% figAxes.Title.String = UIAxes.Title.String;
figAxes.XLabel.String = UIAxes.XLabel.String;
% figAxes.XLim = UIAxes.XLim;
if numel(UIAxes.YAxis) == 1
    ui_children = get(UIAxes,'children');
    copyobj(ui_children, figAxes)
    figAxes.XLim = UIAxes.XLim;
    figAxes.YLim = UIAxes.YLim;
    figAxes.YLabel.String = UIAxes.YLabel.String;
else % numel(UIAxes.YAxis) == 2
    yyaxis(UIAxes,'left')
        ui_children = get(UIAxes,'children');
        lgndNamel = get(ui_children,'displayname');
        yl_lab = UIAxes.YLabel.String;
        yl_lim = UIAxes.YLim;
    yyaxis(figAxes,'left')
        copyobj(ui_children, figAxes)
        figAxes.YLabel.String = yl_lab;
        figAxes.YLim = yl_lim;
    yyaxis(UIAxes,'right')
        ui_children = get(UIAxes,'children');
        lgndNamer = get(ui_children,'displayname');
        yr_lab = UIAxes.YLabel.String;
%         yr_lim = UIAxes.YLim;
    yyaxis(figAxes,'right')
        copyobj(ui_children, figAxes)
        figAxes.YLabel.String = yr_lab;
%         figAxes.YLim = yr_lim;
    lgd = legend({lgndNamel,lgndNamer});
    lgd.Box = UIAxes.Legend.Box;
    lgd.Location = UIAxes.Legend.Location;
end
set(groot, 'Default', struct())
end

