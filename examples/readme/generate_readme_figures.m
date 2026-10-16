function generate_readme_figures
% Reproduce the README figures with the unchanged harmonic package functions.
% Run from the repository root:
%   addpath(fullfile(pwd,'examples','readme')); generate_readme_figures
% All illustrations are synthetic. The time-series example is the first
% stored trial at innovation SD 1, selected before inspecting its projection.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
outdir = fullfile(root,'assets','readme');
builddir = fullfile(root,'build','readme');
mkdir(outdir); mkdir(builddir);
rng(20261016,'twister');
navy = [0.10 0.18 0.29];
teal = [0.00 0.48 0.52];
amber = [0.79 0.38 0.10];
purple = [0.45 0.31 0.65];

% 1. The six-node flow used to define the Live Script's weak-circulation case.
Nodes = [1 0; 0 1; -1 0; 0 -1; 0 2.05; 0 -2.05];
Edges = [1 2; 2 3; 3 4; 4 1; 2 5; 6 4];
X = [2;3;3;2;5;5];
C = flow2matrix(Nodes,Edges,X);
S = adj2simplex(C,1);
B = PH_boundary(S);
[Yg,Yc,Yh] = Hodge_decompose_edgeflow([Edges X],B);
checks.graph.reconstruction = norm(X-Yg-Yc-Yh);
checks.graph.gradient_harmonic_inner_product = abs(Yg'*Yh);
H = flow2matrix(Nodes,Edges,Yh);
hRef = H(sub2ind(size(H),S{2}(:,1),S{2}(:,2)));
checks.graph.harmonic_divergence = norm(B{1}*hRef);
checks.graph.harmonic_expected_error = norm(Yh-[2.5;2.5;2.5;2.5;0;0]);
f = canvas([1440 760]);
headline(f,'Isolate circulation from directional propagation', ...
    'An exact decomposition of the six-node network used in the simulation');
values = {X,Yg,Yh};
names = {'A  Input edge flow','B  Gradient component','C  Harmonic component'};
notes = {'Cycle and two attached branches','Variation explained by node potentials','Persistent circulation around the loop'};
cols = {navy,amber,teal};
for k=1:3
    ax=axes(f,'Position',[0.035+(k-1)*0.33 0.18 0.29 0.61]);
    flowplot(ax,Nodes,Edges,values{k},5,cols{k},true,zeros(0,3));
    title(ax,names{k},'Color',navy,'FontSize',19,'FontWeight','bold');
    label(f,[0.025+(k-1)*0.33 0.105 0.31 0.055],notes{k},14,navy,'center');
end
footer(f,'Input = gradient + harmonic. No filled triangles: curl = 0. Arrow widths share one scale; labels give edge-flow magnitudes.');
print(f,fullfile(outdir,'harmonic-decomposition.png'),'-dpng','-r170');
close(f);

% 2. A full edge-set workflow from the bundled time series, without a
% ground-truth edge mask. The harmonicflow() wrapper accepts lexicographic pairs.
D=load(fullfile(root,'simulationZ_weakcirculation.mat'));
a=find(D.noise_levels==1,1);
Z=D.Zall(:,:,a,1);
TL=TS2timelagged(Z.',1,size(Z,2)-1,1);
A=adj_AntiSymmetric(TL);
Eall=nchoosek(1:size(Z,1),2);
x=A(sub2ind(size(A),Eall(:,1),Eall(:,2)));
result=harmonicflow(x);
Sall=adj2simplex(A,1); Ball=PH_boundary(Sall);
fall=adj2edgeflow(Sall{2}(:,1:2),A);
[gall,call,hall]=Hodge_decompose_edgeflow(fall,Ball);
Hall=flow2matrix((1:size(Z,1))',Sall{2}(:,1:2),hall);
checks.timeseries.wrapper_agreement=norm(double(result.connectivity)-Hall,'fro');
checks.timeseries.harmonic_divergence=norm(Ball{1}*hall);
checks.timeseries.reconstruction=norm(fall(:,3)-gall-call-hall);
checks.timeseries.trial=1;
checks.timeseries.innovation_sd=D.noise_levels(a);
checks.timeseries.size=size(Z);
checks.timeseries.edge_count=size(Eall,1);
f=canvas([1440 730]);
headline(f,'From multivariate time series to a cyclic network component', ...
    'Bundled synthetic data: 6 nodes, 200 samples, first trial at innovation SD = 1');
ax=axes(f,'Position',[0.065 0.27 0.25 0.48]);
hold(ax,'on'); colors=[0.0 0.48 0.52;0.79 0.38 0.10;0.45 0.31 0.65;0.15 0.36 0.66;0.27 0.49 0.25;0.69 0.24 0.31]; ZZ=zscore(Z.');
for k=1:6
    plot(ax,1:200,0.22*ZZ(:,k)+(7-k),'Color',colors(k,:),'LineWidth',1.1);
end
set(ax,'YTick',1:6,'YTickLabel',{'6','5','4','3','2','1'},'FontSize',12, ...
    'Box','off','TickDir','out','XColor',navy,'YColor',navy);
xlim(ax,[1 200]); ylim(ax,[0.25 6.8]);
xlabel(ax,'Time sample'); ylabel(ax,'Node (standardized traces, offset)');
title(ax,'A  Node time series','FontSize',18,'Color',navy);
ax=axes(f,'Position',[0.395 0.27 0.25 0.48]);
imagesc(ax,A); axis(ax,'image'); clim(ax,[-max(abs(A(:))) max(abs(A(:)))]);
colormap(ax,diverging());
set(ax,'XTick',1:6,'YTick',1:6,'FontSize',12,'TickLength',[0 0]);
xlabel(ax,'Target node'); ylabel(ax,'Source node');
title(ax,'B  Directed edge-flow matrix','FontSize',18,'Color',navy);
cb=colorbar(ax,'southoutside'); cb.Label.String='Dominant lagged correlation (signed)'; cb.FontSize=11; cb.Color=navy;
ax=axes(f,'Position',[0.71 0.23 0.26 0.54]);
theta=(0:5)'*2*pi/6+pi/2;
coords=[cos(theta) sin(theta)];
flowplot(ax,coords,Eall,double(result.flows),max(abs(double(result.flows))),teal,false,zeros(0,3));
title(ax,'C  Harmonic projection','FontSize',18,'Color',navy);
label(f,[0.04 0.13 0.29 0.06],'TS2timelagged: correlations at lag 1',13,navy,'center');
label(f,[0.365 0.13 0.30 0.06],'adj_AntiSymmetric: one direction per pair',13,navy,'center');
label(f,[0.69 0.13 0.29 0.06],'harmonicflow: all 15 node pairs are supplied',13,navy,'center');
footer(f,'A workflow example, not a recovery benchmark. Line width encodes absolute harmonic flow; arrowheads show its sign.');
print(f,fullfile(outdir,'timeseries-to-harmonic.png'),'-dpng','-r170');
close(f);

% 3. A filled triangle and an unfilled four-edge loop demonstrate the
% distinction between curl and harmonic flow on a chosen simplicial complex.
N2=[-1.6 0;-1.05 1;-0.5 0;0.6 1;1.7 0;0.6 -1];
E2=[1 2;2 3;1 3;3 4;4 5;5 6;3 6];
reference=[1;1;-1;1;1;1;-1];
C2=flow2matrix(N2,E2,reference);
S2=adj2simplex(C2,2); B2=PH_boundary(S2);
e2=S2{2}(:,1:2);
potential=[-0.25;0.35;0.05;-0.30;0.40;0.20];
G2=B2{1}'*potential;
curl2=1.4*B2{2}(:,1);
cycleMatrix=flow2matrix(N2,E2,[0;0;0;1;1;1;-1]);
harm2=cycleMatrix(sub2ind(size(cycleMatrix),e2(:,1),e2(:,2)));
x2=G2+curl2+harm2;
[g2,c2,h2]=Hodge_decompose_edgeflow([e2 x2],B2);
checks.complex.boundary_chain=norm(B2{1}*B2{2},'fro');
checks.complex.reconstruction=norm(x2-g2-c2-h2);
checks.complex.known_components=max([norm(g2-G2),norm(c2-curl2),norm(h2-harm2)]);
checks.complex.harmonic_divergence=norm(B2{1}*h2);
checks.complex.harmonic_curl=norm(B2{2}'*h2);
checks.complex.orthogonality=max(abs([g2'*c2,g2'*h2,c2'*h2]));
f=canvas([1440 1020]);
headline(f,'Distinguish local triangular circulation from an unfilled loop', ...
    'The chosen simplicial complex determines which cycles are curl and which are harmonic');
v2={x2,g2,c2,h2}; col2={navy,amber,purple,teal};
titles={'A  Input flow','B  Gradient: node-potential variation', ...
    'C  Curl: circulation around a filled face','D  Harmonic: circulation around the open loop'};
for k=1:4
    row=ceil(k/2); col=mod(k-1,2);
    ax=axes(f,'Position',[0.04+col*0.50 0.51-(row-1)*0.39 0.42 0.30]);
    flowplot(ax,N2,e2,v2{k},max(abs(x2)),col2{k},false,S2{3}(:,1:3));
    title(ax,titles{k},'Color',navy,'FontSize',17,'FontWeight','bold');
end
footer(f,'Input = gradient + curl + harmonic. Shading marks the filled triangle. All four panels use the same arrow-width scale.');
print(f,fullfile(outdir,'hodge-components.png'),'-dpng','-r170');
close(f);

% Numerical checks are for these examples, not an end-to-end package test.
assert(checks.graph.reconstruction<1e-10);
assert(checks.graph.harmonic_expected_error<1e-10);
assert(checks.timeseries.wrapper_agreement<1e-6);
assert(checks.timeseries.harmonic_divergence<1e-10);
assert(checks.complex.known_components<1e-10);
assert(checks.complex.orthogonality<1e-10);
checks.matlab=version;
checks.source_commit='04e0574a79c87c94c7c702a5446ab8744a7d2a57';
fid=fopen(fullfile(builddir,'numerical-checks.json'),'w');
fwrite(fid,jsonencode(checks,PrettyPrint=true)); fclose(fid);
save(fullfile(builddir,'figure-data.mat'),'X','Yg','Yc','Yh','Nodes','Edges', ...
    'Z','A','result','N2','e2','x2','g2','c2','h2','checks');
disp(checks);
fprintf('README figures exported successfully.\n');
end

function f=canvas(sz)
f=figure('Visible','off','Color','white','Position',[80 80 sz]);
set(f,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesColor','white','DefaultAxesXColor',[0.10 0.18 0.29], ...
    'DefaultAxesYColor',[0.10 0.18 0.29],'PaperPositionMode','auto','InvertHardcopy','off');
end

function headline(f,heading,subheading)
label(f,[0.04 0.925 0.92 0.055],heading,24,[0.10 0.18 0.29],'left');
label(f,[0.04 0.875 0.92 0.04],subheading,15,[0.35 0.41 0.49],'left');
end

function footer(f,str)
annotation(f,'line',[0.04 0.96],[0.078 0.078],'Color',[0.84 0.87 0.90]);
label(f,[0.04 0.016 0.92 0.049],str,13,[0.35 0.41 0.49],'left');
end

function label(f,pos,str,sz,color,align)
annotation(f,'textbox',pos,'String',str,'EdgeColor','none','FontName','Arial', ...
    'FontSize',sz,'Color',color,'HorizontalAlignment',align,'VerticalAlignment','middle','Interpreter','none');
end

function flowplot(ax,nodes,edges,weights,scale,color,showWeights,faces)
hold(ax,'on');
for k=1:size(faces,1)
    patch(ax,nodes(faces(k,:),1),nodes(faces(k,:),2),[0.91 0.88 0.96],'EdgeColor','none');
end
for k=1:size(edges,1)
    xy=nodes(edges(k,:),:);
    plot(ax,xy(:,1),xy(:,2),'-','Color',[0.87 0.89 0.91],'LineWidth',1.1);
end
for k=1:size(edges,1)
    w=weights(k);
    if abs(w)<1e-9, continue; end
    ij=edges(k,:); if w<0, ij=fliplr(ij); end
    a=nodes(ij(1),:); b=nodes(ij(2),:);
    u=(b-a)/norm(b-a); v=[-u(2) u(1)];
    a=a+0.10*u; b=b-0.11*u;
    lw=0.65+4.3*abs(w)/scale;
    plot(ax,[a(1) b(1)],[a(2) b(2)],'-','Color',color,'LineWidth',lw);
    head=0.11+0.035*abs(w)/scale;
    tip=b; base=b-head*u;
    P=[tip;base+0.46*head*v;base-0.46*head*v];
    patch(ax,P(:,1),P(:,2),color,'EdgeColor',color);
    if showWeights
        mid=(a+b)/2+0.12*v;
        text(ax,mid(1),mid(2),sprintf('%.1f',abs(w)),'FontSize',13,'FontWeight','bold', ...
            'Color',color,'HorizontalAlignment','center','BackgroundColor','white','Margin',1);
    end
end
scatter(ax,nodes(:,1),nodes(:,2),330,'white','filled','MarkerEdgeColor',[0.10 0.18 0.29],'LineWidth',1.5);
for k=1:size(nodes,1)
    text(ax,nodes(k,1),nodes(k,2),num2str(k),'FontSize',13,'FontWeight','bold', ...
        'Color',[0.10 0.18 0.29],'HorizontalAlignment','center','VerticalAlignment','middle');
end
axis(ax,'equal'); axis(ax,'off');
xlim(ax,[min(nodes(:,1))-0.35 max(nodes(:,1))+0.35]);
ylim(ax,[min(nodes(:,2))-0.30 max(nodes(:,2))+0.30]);
end

function c=diverging
left=[0.22 0.35 0.67]; mid=[0.98 0.98 0.98]; right=[0.77 0.25 0.17];
t=linspace(0,1,128)'; c=[left.*(1-t)+mid.*t;mid.*(1-t)+right.*t];
end
