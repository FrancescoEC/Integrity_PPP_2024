function  []= saveallfigures(tempdir)
if nargin==1
FolderName = tempdir;   % Your destination folder
else
    FolderName = pwd;
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
  saveas(FigHandle,fullfile(FolderName, [FigName '.png']));
end
end
