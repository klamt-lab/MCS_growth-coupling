function p = plot_mcs_relationships(reac_names,mcs_weakGC,mcs_strongGC,mcs_substUp)

[~,~,~,a] = compare_mcs_sets(mcs_strongGC,mcs_weakGC);
[~,~,~,b] = compare_mcs_sets(mcs_substUp,mcs_weakGC);
[~,~,~,c] = compare_mcs_sets(mcs_substUp,mcs_strongGC);
e = ~(c*a) .* b;

% check for completeness
if all([max(a,[],2)>=2; 1]) && all([max(b,[],2)>=2; 1])
    disp('Hierarchy of MCS verified.');
% check for correctness
elseif all([max(a,[],2)~=1; 1]) && all([max(b,[],2)~=1; 1]) && all([max(c,[],2)~=1; 1])
    disp('Hierarchy of MCS verified, but some MCS must not have been found.');
else
    disp('No hierarchy of MCS.');
end

link_weak_strong            = cell(size(mcs_weakGC,2),1);
link_weak_strong_type       = zeros(size(mcs_weakGC,2),1);
[ro,co] = find(a);
for i = unique(co(:))'
    link_weak_strong{i} = ro(co==i)';
    link_weak_strong_type(i) = max(a(:,i));
end
link_strong_substrate       = cell(size(mcs_strongGC,2),1);
link_strong_substrate_type  = zeros(size(mcs_strongGC,2),1);
[ro,co] = find(c);
for i = unique(co(:))'
    link_strong_substrate{i} = ro(co==i)';
    link_strong_substrate_type(i) = max(c(:,i));
end
link_weak_substrate         = cell(size(mcs_weakGC,2),1);
link_weak_substrate_type    = zeros(size(mcs_weakGC,2),1);
[ro,co] = find(e);
for i = unique(co(:))'
    link_weak_substrate{i} = ro(co==i)';
    link_weak_substrate_type(i) = max(e(:,i));
end

% draw graph
G = digraph();
x = [];
y = [];
mcs_text = cell(size(mcs_weakGC,2)+size(mcs_strongGC,2)+size(mcs_substUp,2),1);
for i = 1:size(mcs_weakGC,2)
    G = G.addnode(['weak_GC' num2str(i)]);
%     mcs_text{G.numnodes,1} = num2str(sum(abs(mcs_weakGC(:,i))));
%     mcs_text{G.numnodes,1} = mcs2text(reac_names,mcs_weakGC(:,i));
    mcs_text{G.numnodes,1} = [num2str(sum(abs(mcs_weakGC(:,i)))) ' (' mcs2text(reac_names,mcs_weakGC(:,i)) ')'];
    x = [x 0];
    y = [y size(mcs_weakGC,2)-i];
end
for i = 1:size(mcs_strongGC,2)
    G = G.addnode(['strong_GC' num2str(i)]);
%     mcs_text{G.numnodes,1} = num2str(sum(abs(mcs_strongGC(:,i))));
%     mcs_text{G.numnodes,1} = mcs2text(reac_names,mcs_strongGC(:,i));
    mcs_text{G.numnodes,1} = [num2str(sum(abs(mcs_strongGC(:,i)))) ' (' mcs2text(reac_names,mcs_strongGC(:,i)) ')'];
%     x = [x 3+i];
%     y = [y size(mcs_strongGC,2)-i];
    x = [x 5];
    y = [y max(size(mcs_weakGC,2),size(mcs_substUp,2))+size(mcs_weakGC,2)-i];
end
for i = 1:size(mcs_substUp,2)
    G = G.addnode(['substrate_uptake_C' num2str(i)]);
%     mcs_text{G.numnodes,1} = num2str(sum(abs(mcs_substUp(:,i))));
%     mcs_text{G.numnodes,1} = mcs2text(reac_names,mcs_substUp(:,i));
    mcs_text{G.numnodes,1} = [num2str(sum(abs(mcs_substUp(:,i)))) ' (' mcs2text(reac_names,mcs_substUp(:,i)) ')'];
%     x = [x 3+i];
%     y = [y -3];
    x = [x 10];
    y = [y size(mcs_weakGC,2)-i];
end
mcs_text = strrep(mcs_text,'EX\_o2\_e','O_2');

edges = nan(0,3); % source, target, type (identical = 3, subset = 2, superset = 1)
for i = 1:size(mcs_weakGC,2)
    target = link_weak_strong{i};
    for j = target
        edges = [edges; i j+size(mcs_weakGC,2) link_weak_strong_type(i)];
    end
    target = link_weak_substrate{i};
    for j = target
        edges = [edges; i j+size(mcs_weakGC,2)+size(mcs_strongGC,2) link_weak_substrate_type(i)];
    end
end
for i = 1:size(mcs_strongGC,2)
    target = link_strong_substrate{i};
    for j = target
        edges = [edges; i+size(mcs_weakGC,2) j+size(mcs_weakGC,2)+size(mcs_strongGC,2) link_strong_substrate_type(i)];
    end
end
color = nan(size(edges,1),1);
color(edges(:,3) == 1) = 8;
color(edges(:,3) == 2) = 1.5;
color(edges(:,3) == 3) = 0;
for i = 1:size(edges,1)
    if edges(:,3) ~= 1
        G = G.addedge(edges(i,1),edges(i,2));
    else
        G = G.addedge(edges(i,2),edges(i,1));
    end
end

% plot
p = plot(G,'XData',x,'YData',y,'EdgeCData',color,'NodeLabel',mcs_text);

% title
plot_title = {};
if ~isempty(mcs_weakGC)
    plot_title = [plot_title {'weak GC MCS'}];
end
if ~isempty(mcs_strongGC)
    plot_title = [plot_title {'strong GC MCS'}];
end
if ~isempty(mcs_substUp)
    plot_title = [plot_title {'subup C MCS'}];
end
title(strjoin(plot_title,'     vs     '));

% colors
% colormap('Turbo');
caxis([0 10]);

% layout
if ~isempty(mcs_weakGC)
    sources = 1:size(mcs_weakGC,2);
elseif ~isempty(mcs_strongGC) && ~isempty(mcs_substUp)
    sources = 1:size(mcs_strongGC,2);
else
    sources = [];
end
if ~isempty(mcs_substUp)
    sinks = size(mcs_weakGC,2)+size(mcs_strongGC,2)+(1:size(mcs_substUp,2));
elseif ~isempty(mcs_strongGC) && ~isempty(mcs_weakGC)
    sinks = size(mcs_weakGC,2)+(1:size(mcs_strongGC,2));
else
    sinks = [];
end
if ~isempty(sinks) && ~isempty(sources)
    layout(p,'layered','Direction','right','Sources',sources,'Sinks',sinks);
elseif isempty(sinks) && ~isempty(sources)
    layout(p,'layered','Direction','right','Sources',sources);
elseif ~isempty(sinks) && isempty(sources)
    layout(p,'layered','Direction','right','Sinks',sinks);
else
    layout(p,'layered','Direction','right');
end
% Prepare output and colors
p.NodeColor = 'black';
p.EdgeColor = 'black';
set(gca,'visible','off');

end

function ko_ki_text = mcs2text(reac_names,mcs)
    reac_names = strrep(cellstr(reac_names),'_','\_');
    kos = find(mcs<0);
    kis = find(mcs>0);
    if ~isempty(kis) 
        kis = [strjoin(strcat('+',reac_names(kis)),', ') ', '];
    else
        kis = [];
    end
    if ~isempty(kos) 
        kos = strjoin(reac_names(kos,:),', ');
    else
        kos = [];
    end
    ko_ki_text = [kis kos];
end
