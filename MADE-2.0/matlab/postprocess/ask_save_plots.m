function files = ask_save_plots(figs, default_folder)
%ASK_SAVE_PLOTS Ask whether to save the supplied figures; keep them open.
%   ask_save_plots() prompts for all currently open figures.
%   ask_save_plots(mech_fig) prompts for the surrogate FEM figure only.
%   Outputs: PNG at 300 dpi and an editable FIG (OFIG with Octave fallback).
%   No solver is rerun. Existing files are never overwritten.
if nargin<1, figs=flipud(findall(0,'Type','figure')); end
if nargin<2 || isempty(default_folder)
    default_folder=fullfile(pwd,['plots_' datestr(now,'yyyymmdd_HHMMSS')]);
end
files={};
% Figures may have been closed by the user since they were created.
keep=false(size(figs));
for k=1:numel(figs)
    keep(k)=ishghandle(figs(k)) && strcmp(get(figs(k),'Type'),'figure');
end
figs=figs(keep);
if isempty(figs), fprintf('No open plots to save.\n'); return; end
answer=strtrim(input(sprintf('Save the %d open plot figure(s), including surrogate FEM if run? [y/N]: ',numel(figs)),'s'));
if ~strcmpi(answer,'y'), return; end
folder=strtrim(input(sprintf('Output folder [%s]: ',default_folder),'s'));
if isempty(folder),folder=default_folder;end
if ~exist(folder,'dir')
    [ok,msg]=mkdir(folder);
    if ~ok, warning('ask_save_plots:folder','Could not create output folder: %s',msg);return;end
end
for k=1:numel(figs)
    f=figs(k);
    if ~ishghandle(f),continue;end
    name=get(f,'Name');
    if isempty(name),name='plot';end
    name=regexprep(char(name),'[^A-Za-z0-9_-]+','_');
    name=name(1:min(numel(name),80));
    if isempty(name),name='plot';end
    base=fullfile(folder,sprintf('%03d_%s',k,name));
    candidate=base;suffix=0;
    while exist([candidate '.png'],'file') || exist([candidate '.fig'],'file') || exist([candidate '.ofig'],'file')
        suffix=suffix+1;candidate=sprintf('%s_%03d',base,suffix);
    end
    % Save each format independently so a graphics-driver failure does not
    % prevent saving the editable figure or any subsequent plot.
    try
        png=[candidate '.png'];print(f,png,'-dpng','-r300');
        files{end+1}=png;fprintf('Saved: %s\n',png); %#ok<AGROW>
    catch err
        warning('ask_save_plots:png','PNG export failed for figure %d: %s',k,err.message);
    end
    try
        if exist('savefig','file') || exist('savefig','builtin')
            editable=[candidate '.fig'];savefig(f,editable);
        else
            editable=[candidate '.ofig'];hgsave(f,editable);
        end
        files{end+1}=editable;fprintf('Saved: %s\n',editable); %#ok<AGROW>
    catch err
        warning('ask_save_plots:editable','Editable export failed for figure %d: %s',k,err.message);
    end
end
fprintf('%d plot file(s) saved. Figures remain open.\n',numel(files));
end
