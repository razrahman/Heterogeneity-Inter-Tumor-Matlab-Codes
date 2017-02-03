function [vLineHandleTree hLineHandleTree textHandleTree] = plot(obj, heightTree, varargin)
%% PLOT  Plot the tree.
% 
%   PLOT(T) lay out the tree T on new axes, compying to Edward Tufte
%   recommandations.
%
%   PLOT(T, LT), where LT is a synchronized tree made of scalars, lay out
%   the tree, using LT data to specify the length of each vertical branch.
%   Use the empty array [] to use a default of 1 for all branches.
%
%   PLOT(T, LT, 'PropertyName', PropertyValue, ...) allows to specify extra
%   parameters for the plot:
%
%       'Ylabel' - a string or a cell array of strings: Use a label for the
%       Y axis. The axis itself becomes visible and ticks are drawn on the
%       tree.
%
%       'X' - a scalar: X leftmost position of the tree.
%
%       'Width' - a scalar: Width of the tree.
%
%       'TextRotation' - a scalar: Rotation, in degrees, of the label that
%       is printed above each node.
%
%       'Parent - an axes handle: The axes handle to plot the tree in. The
%       axes are not changed, except for the YLabel.
%
%   [ VL HL TH ] = PLOT(T, ...) returns three synchronized trees containing
%   respectively the handles for the vertical lines, the horizontal lines
%   and the text labels of each node.
%
%   EXAMPLE:
%   [ lineage duration ] = tree.example; % 1st one is made of strings only, 2nd one of integers
%   slin = lineage.subtree(19); % Work on a subset
%   sdur = duration.subtree(19);
%   [vlh hlh tlh] = slin.plot(sdur, 'YLabel', {'Division time' '(min)'});
%   rcolor = [ 0.6 0.2 0.2 ];
%   aboveTreshold = sdur > 10; % true if longer than 10 minutes
%   iterator = aboveTreshold.depthfirstiterator;
%   for i = iterator
%    if  aboveTreshold.get(i)
%        set( vlh.get(i), 'Color' , rcolor )
%        set( hlh.get(i), 'Color' , rcolor )
%        set( tlh.get(i), 'Color' , rcolor )
%    end
% end

% Jean-Yves Tinevez <tinevez AT pasteur DOT fr> March 2012

    %% CONSTANTS
    
    LINE_COLOR = [ 0.3 0.3 0.3 ];
    
    %% Deal with input
    
    if nargin < 2 || isempty(heightTree)
        heightTree = tree(obj, 1);
    end
    
    parser = inputParser;
    parser.addParamValue('YLabel', [], @(x) ischar(x) || iscell(x));
    parser.addParamValue('TextRotation', 0, @(x) isnumeric(x) && isscalar(x) );
    parser.addParamValue('Parent', [], @ishandle);
    parser.addParamValue('X', 0, @(x) isnumeric(x) && isscalar(x) );
    parser.addParamValue('Width', 1, @(x) isnumeric(x) && isscalar(x) );
    
    parser.parse(varargin{:});
    ylbl    = parser.Results.YLabel;
    textrot = mod(parser.Results.TextRotation, 360);
    ax      = parser.Results.Parent;
    xcorner = parser.Results.X;
    xwidth  = parser.Results.Width;
    
    %% Compute the column width
    
    width = tree(obj, 'clear');

    % Put 1 at the leaves
    iterator = width.depthfirstiterator;
    for i = iterator
       if width.isleaf(i)
           width = width.set(i, 1);
       end
    end
    
    % Cumsum
    width = width.recursivecumfun(@sum);
    
    % Normalize
    maxWidth = width.get(1);
    width = width .* ( xwidth / maxWidth );
    
    %% Compute the X *column* width
    
    xcol = tree(width, 'clean');
    xcol = xcol.set(1, 0);
    
    previous = 0;
    parent = 1;
    iterator = width.breadthfirstiterator;
    
    for i = iterator(2 : end) % The root is already done
       
        newParent = xcol.getparent(i);
        if newParent ~= parent
           % We just changed branch
           parent = newParent;
           previous = xcol.get(parent);
        end
        
        w = width.get(i);
        xcol = xcol.set(i, previous);
        
        previous = previous + w;
    end
    
    %% Compute the actual X position
  
    xpos = tree(width, 'clean');
    xpos = xpos.set(1, 1);
    iterator = xpos.breadthfirstiterator;
    for i = iterator
        xpos = xpos.set(i, xcol.get(i) + width.get(i)/2);
    end
    
    % Max of x position
    maxXpos = -1;
    for i = iterator
        xp = xpos.get(i);
        if xp > maxXpos
            maxXpos = xp;
        end
    end
    
    
    
    %% Compute the Y position
    
    ypos = tree(obj, 'clear');
    ypos = ypos.set(1, heightTree.get(1));
    iterator = ypos.depthfirstiterator;
    iterator(1) = []; % Skip the root
    
    maxHeight = heightTree.get(1);
    
    for i = iterator
       parent = ypos.getparent(i);
       parentPos = ypos.get(parent);
       height = heightTree.get(i);
       ypos = ypos.set(i, parentPos + height);
       
       if maxHeight < parentPos + height
           maxHeight = parentPos + height;
       end
    end
    
    %% Prepare the axes
    
    if isempty(ax) 
        ax = axes( ...
            'FontName', 'Courier new', ...
            'FontSize', 9, ...
            'YDir', 'reverse', ...
            'TickDir', 'out', ...
            'XTickLabel', '', ...
            'XTick', [], ...
            'XLim', [0 maxXpos * 1.05], ...
            'XColor', 'w');
    end
    
    if isempty(ylbl)
        set(ax, ...
            'YColor', 'w', ...
            'YTick', [], ...
            'YTickLabel', '')
    else
        ylabel(ylbl, ...
            'HorizontalAlignment', 'right', ...
            'Rotation', 0)
    end
    hold(ax, 'on')
    
    %% A first iteration for the vertical bars
        
    % Prepare holder for the vertical line handles
    vLineHandleTree = tree(obj, 'clear');
    
    iterator = obj.depthfirstiterator;
    for i = iterator
        
        % Vertical bars -> to parent
        
        y1 = ypos.get(i);
        
        if isempty(y1)
            continue
        end
        
        y2 = y1 - heightTree.get(i);
        
        x1 = xpos.get(i) + xcorner;
        x2 = x1;
        
        hl = line([x1 x2], [y1 y2], ...
            'Color', LINE_COLOR, ...
            'LineWidth', 1);
        
        vLineHandleTree = vLineHandleTree.set(i, hl);
        
    end
    
    % If we were given a height tree, draw white ticks on the tree, a la
    % Tufte.

    if nargin >= 2
        xl = xlim;
        ticks = get(ax, 'YTick');
        
        if ~isempty(ylbl)
            % Chop to min and max, a la Tufte
            ticks(end) = maxHeight;
            set(ax, 'YTick', ticks, 'YLim', [0 maxHeight])
        end
        
        % White ticks
        bgColor = get(ax, 'Color');
        for t = ticks
            line( xl + xl(2)/1e2, [t t], ...
                'LineWidth', 2, ...
                'Color', bgColor)
        end
    end
    
    %% New iteration for the bars and the content
        
    % Prepare the holder for the text handles
    textHandleTree = tree(obj, 'clear');
    
    % Prepare the holder for horizontal line handles
    hLineHandleTree = tree(obj, 'clear');
    
    % Prepare display of text
    if textrot <  45 || (textrot >  135 && textrot < 225) || textrot > 315
        halign = 'center';
        valign = 'middle';
        contentfun = @(x) { x ' ' ' ' };
    else
        halign = 'left';
        valign = 'middle';
        contentfun = @(x) [ ' ' x ];
    end
    
    for i = iterator
        
         y1 = ypos.get(i);
        
        if isempty(y1)
            continue
        end
        
        y2 = y1 - heightTree.get(i);
        
        x1 = xpos.get(i) + xcorner;
        
        % The label = content
        content = obj.get(i);
        if isempty(content)
            content = '�';
        end
        if ~ischar(content)
            content = num2str(content);
        end
        
        ht = text(x1, y2, contentfun(content), ...  A hack to have text displayed above bars
            'HorizontalAlignment', halign,...
            'Rotation', textrot, ...
            'VerticalAlignment', valign, ...
            'FontName', 'Courier new', ...
            'Interpreter', 'none', ...
            'FontSize', 12);
        
        textHandleTree = textHandleTree.set(i, ht);
        
        % Horizontal bars -> children
        if obj.isleaf(i)
            continue
        end
        
        children = obj.getchildren(i);
        allX = cell2mat(xpos.Node(children)) + xcorner;
        
        y2 = y1;
        x1 = min(allX);
        x2 = max(allX);

        hl = line([x1 x2], [y1 y2], ...
            'Color', LINE_COLOR, ...
            'LineWidth', 5);
        
        hLineHandleTree = hLineHandleTree.set(i, hl);
        
    end
    
end