function h=myPlot(x,y,symb,Color,xlb,ylb,legnd,xtk,ytk)
% myPlot(x,y,symb,Color,xlb,ylb,legnd,xtk,ytk)
if ~exist('Color','var') || isempty(Color)
    
    Color = lines(size(y,2));
    
end

if size(Color,1) <size(y,2)
    Color = repmat(Color,[ceil(size(y,2)/size(Color,1)),1]);
end
 symb=str2cell(symb);
if (size(x,1)==1||size(x,2)==1) && (size(y,1)>1&&size(y,2)>1)

  
    if length(symb)==1
        symb=repmat(symb,[1,size(y,2)]);
    end
    for i=1:size(y,2)
        gco= plot(x,y(:,i),symb{i},'LineWidth',1,'Color',Color(i,:));
      %  set(gco,'Color',Color(i,:));
        hold on;
        h(i)=gco;
    end
elseif ~any(size(x)~=size(y)) && (size(y,1)>1&&size(y,2)>1)
    
    if length(symb)==1
        symb=repmat(symb,[1,size(y,2)]);
    end
    
    for i=1:size(y,2)
        gco= plot(x(:,i),y(:,i),symb{i},'LineWidth',1);
        set(gco,'Color',Color(i,:));
        hold on;
        h(i)=gco;
    end
else
   h= plot(x,y,symb{1},'LineWidth',1);
    set(h,'Color',Color(1,:));
end
set(gca,'FontSize',12);
if exist('xlb','var')
    xlabel(xlb);
end

if exist('ylb','var') 
    ylabel(ylb);
end

if exist('legnd','var') && ~isempty(legnd)
    legend(legnd);
end

if exist('xtk','var')  && ~isempty(xtk)
    xlim(min_max(xtk));
    if length(xtk)>2
        set(gca,'XTick',xtk);
    end
    
end

if exist('ytk','var') && ~isempty(ytk)
    ylim(min_max(ytk));
     if length(ytk)>2
         
        set(gca,'YTick',ytk);
    end
end

 set(gca,'TickLength',[0.02,0.02]);
 box on;
