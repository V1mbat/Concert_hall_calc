function [] = exportPlot(f,filePathP,plotName,varargin)
    %Export plots as pdf, 
    %f - Plot variable f = figure(...)
    %filePathP = fullfile(pwd,"Plots/")
    %plotName = 'exp_X.y_NameForYourPlot
    n = nargin;
    lcn = ['north','south','east','west','northeast','northwest','southeast',...
        'southwest','northoutside','southoutside','eastoutside','westoutside',...
        'northeastoutside','northwestoutside','southeastoutside','southwestoutside',...
        'best','bestoutside','layout','none'];
    for i = 1:size(varargin,2)
    if isa(varargin{i},"logical")
    plotLogic = varargin{i};
    elseif isa(varargin{i},"char") && logical(contains(lcn,varargin{i}))
    loc =  varargin{i};
    else 
        loc = 'bestoutside';
    end
    end
    
    %Fix plot sizeand legend
    thickenall_big
    f.Position = [300, 0, 990, 788.8]; 
    if  ~isempty(f.CurrentAxes.Legend) %|| 
    f.CurrentAxes.Legend.Location=loc;        
    end
    f.CurrentAxes.GridLineStyle ="-";
    grid on

    if plotLogic
    %Setup save name
    plotName = strcat(plotName,'.pdf');
    saveNameP = fullfile(filePathP,plotName);
    %Save the plot
    exportgraphics(f,saveNameP,"Resolution",600,"BackgroundColor","none")
    end
end