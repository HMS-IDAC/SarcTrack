classdef timeLapseViewTool < handle    
    properties
        Figure
        FigureContext
        Axis
        PlaneHandle
        PlaneIndex
        NPlanes
        Volume
        Dialog
        LowerThreshold
        UpperThreshold
        LowerThresholdSlider
        UpperThresholdSlider
        FigurePlot
        AxisPlotPosition
        Data
    end
    
    methods
        function tool = timeLapseViewTool(V,varargin)
% timeLapseViewTool(V)
% A tool to visualize a time lapse movie
% V should be m by n by z 'double' and in the range [0,1]
% z is the number of frames; for small z (e.g. 1...5), use stackViewTool instead
%
% example:
% load mri
% V = double(squeeze(D))/255;
% timeLapseViewTool(V)

            tool.Volume = V;
            tool.NPlanes = size(V,3);
            tool.PlaneIndex = 1;
            
            tool.LowerThreshold = 0;
            tool.UpperThreshold = 1;
               
            % z
            tool.Figure = figure('Name','Frame','NumberTitle','off','CloseRequestFcn',@tool.closeTool);
            tool.Axis = axes('Parent',tool.Figure,'Position',[0 0 1 1]);
            
            I = tool.Volume(:,:,tool.PlaneIndex);
            tool.PlaneHandle = imshow(tool.applyThresholds(I));
            
            tool.Axis.Title.String = sprintf('frame %d', tool.PlaneIndex);
            
            dwidth = 1000;%300;
            dborder = 10;
            cwidth = dwidth-2*dborder;
            cheight = 20;
            
            tool.Dialog = dialog('WindowStyle', 'normal', 'Resize', 'on',...
                                'Name', 'TimeLapseViewTool',...
                                'CloseRequestFcn', @tool.closeTool,...
                                'Position',[100 100 dwidth 4*dborder+4*cheight]);
                            
            % z slider
            uicontrol('Parent',tool.Dialog,'Style','text','String','f','Position',[dborder 3*dborder+3*cheight 20 cheight]);
            slider = uicontrol('Parent',tool.Dialog,'Style','slider','Min',1,'Max',tool.NPlanes,'Value',tool.PlaneIndex,'Position',[dborder+20 3*dborder+3*cheight cwidth-20 cheight],'Tag','zs');
            addlistener(slider,'Value','PostSet',@tool.continuousSliderManage);

            % lower threshold slider
            uicontrol('Parent',tool.Dialog,'Style','text','String','_t','Position',[dborder 2*dborder+cheight 20 cheight]);
            tool.LowerThresholdSlider = uicontrol('Parent',tool.Dialog,'Style','slider','Min',0,'Max',1,'Value',tool.LowerThreshold,'Position',[dborder+20 2*dborder+cheight cwidth-20 cheight],'Tag','lts');
            addlistener(tool.LowerThresholdSlider,'Value','PostSet',@tool.continuousSliderManage);
            
            % upper threshold slider
            uicontrol('Parent',tool.Dialog,'Style','text','String','^t','Position',[dborder dborder 20 cheight]);
            tool.UpperThresholdSlider = uicontrol('Parent',tool.Dialog,'Style','slider','Min',0,'Max',1,'Value',tool.UpperThreshold,'Position',[dborder+20 dborder cwidth-20 cheight],'Tag','uts');
            addlistener(tool.UpperThresholdSlider,'Value','PostSet',@tool.continuousSliderManage);
            
            if nargin > 1
                tool.FigurePlot = figure('Name','Frame','NumberTitle','off','CloseRequestFcn',@tool.closeTool);
                tool.Data = varargin{1};
                plot(1:length(tool.Data),tool.Data)
                hold on
                tool.AxisPlotPosition = plot(tool.PlaneIndex,tool.Data(tool.PlaneIndex),'o');
                hold off
            end
            
%             uiwait(tool.Dialog)
        end
        
        function continuousSliderManage(tool,~,callbackdata)
            tag = callbackdata.AffectedObject.Tag;
            value = callbackdata.AffectedObject.Value;
            if strcmp(tag,'uts') || strcmp(tag,'lts')
                if strcmp(tag,'uts')
                    tool.UpperThreshold = value;
                elseif strcmp(tag,'lts')
                    tool.LowerThreshold = value;
                end
                
                I = tool.Volume(:,:,tool.PlaneIndex);
                tool.PlaneHandle.CData = tool.applyThresholds(I);
            elseif strcmp(tag,'zs')
                tool.PlaneIndex = round(value);

                I = tool.Volume(:,:,tool.PlaneIndex);
                tool.PlaneHandle.CData = tool.applyThresholds(I);

                tool.Axis.Title.String = sprintf('frame %d', tool.PlaneIndex);
                if ~isempty(tool.FigurePlot)
                    tool.AxisPlotPosition.XData = tool.PlaneIndex;
                    tool.AxisPlotPosition.YData = tool.Data(tool.PlaneIndex);
                end
            end
        end
        
        function T = applyThresholds(tool,I)
            T = I;
            T(T < tool.LowerThreshold) = tool.LowerThreshold;
            T(T > tool.UpperThreshold) = tool.UpperThreshold;
            T = T-min(T(:));
            T = T/max(T(:));
        end
        
        function closeTool(tool,~,~)
            delete(tool.Figure)
            delete(tool.FigurePlot);
            delete(tool.Dialog);
        end
    end
end
