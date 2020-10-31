size_pic=11;
pic = findobj(gcf); 
allaxes  = findall(pic,'Type','axes');
alltext  = findall(pic,'Type','text');
set(allaxes,'FontSize', size_pic, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal');
set(alltext,'FontSize', size_pic, 'FontName', 'Times New Roman', ...
            'FontWeight', 'normal', 'FontAngle', 'normal','Interpreter', 'latex');
set(gca,'TickDir','in'); 
