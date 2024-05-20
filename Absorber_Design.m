
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     Author: Abhishek Gupta                          %%%
%%%                    Title: Research Assistant                        %%%
%%%                    email: agupta253@wisc.edu                        %%%
%%%                Advisor: Ramathasan Thevamaran                       %%%
%%%               Department of Mechanical Engineering                  %%%
%%%                University of Wisconsin-Madison                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
set(0,'defaultfigurewindowstyle','docked')
addpath('F:\MATLAB Tools\Colormaps\Colormaps (5)\Colormaps')
set(groot,'defaultAxesFontSize',18)
set(0,'defaultfigurecolor',[1 1 1])
set(0, 'DefaultLineLineWidth', 3);
warning off
Colors = ([0.5, 0, 0.5;0, 0.5, 0.5;0.5, 0.5, 0;1, 0.5, 0.5;0.529, 0.808, 0.922;1, 0.5, 0.314; 0.902, 0.902, 0.18; 1,0.2,0.3; 0,0.2,0.8]);

%% Input
%{
mas = input(['Enter the mass of accelerating object (m) in kg','\n']);
ac = input(['Enter the peak acceleration (ac) in m/s^2','\n']);
vi = input(['Enter initial velocity (vi) in m/s','\n']);
Are_t = input(['Enter a targeted cross section area in m^2','\n']);
%}
mas = 2;
ac = 1000;
vi = 4;
Are_t = 0.01;
%}
%%
Msth = readtable('samp_dat.xlsx');      % foam material data
Msth2 = table2array(Msth(:,[2:10,13:16])); 
sz = size(Msth2);

for ind =1:sz(1)

    filename = [char(table2cell(Msth(ind,12))),'.mat'];
    load(filename);
    suf = [char(table2cell(Msth(ind,11)))];
    sampn = Msth2(ind,11);
    
    Strain_cell = evalin('base', [suf,'_Avg_Strain_L']);
    Stress_cell = evalin('base', [suf,'_Avg_Stress_L']);

    Strain_L{ind} = Strain_cell{sampn};          % Loading strain
    Stress_L{ind} = Stress_cell{sampn};          % Loading stress in MPa

    smoothness_parameter = 0.99999;  
    fitting_spline = csaps(Strain_L{ind}, Stress_L{ind}, smoothness_parameter);
    intervals{ind} = linspace(min(Strain_L{ind}), max(Strain_L{ind}), 40000);  
    fitting_points{ind} = fnval(fitting_spline, intervals{ind});

    Epd2(ind) = Msth2(ind,13);           % Critical strain
    %Epd2(ind) = 0.5;

    figure(ind)
    hold on
    plot(Strain_L{ind},Stress_L{ind})
    plot(intervals{ind},fitting_points{ind},'LineStyle','--')
    
    dens(ind) = 1000*Msth2(ind,5);  % density in kg/m^3
    rhob(ind) = Msth2(ind,5)/Msth2(ind,10);  % relative density

    % Fitting power law
    cos = Epd2(ind); % max Strain for fitting
    [mn,In] = min(abs(intervals{ind}-cos));
    xfit = intervals{ind}(100:In);
    yfit = fitting_points{ind}(100:In);
    f1 = fit(xfit.',yfit.','El*x^ll');
    Elv(ind) = 1e6*f1.El;  % fit modulus in pa = El*(lam+1)^beta
    lam(ind) = f1.ll;
    
    figure(ind)
    xplt = linspace(0,Epd2(ind),1000);
    plot(xplt,1e-6*Elv(ind)*(xplt.^lam(ind)),'LineStyle','--')
    box on
    xlabel('Strain','Interpreter','Latex','Rotation',0,'FontName','Arial','color','k','FontSize',24)
    ylabel('Stress (MPa)','Interpreter','Latex','Rotation',90,'FontName','Arial','color','k','FontSize',24)
    set(gca,'Xscale','linear','LineWidth',1.5)
    xlim([0,1])
    %xline(Epd2(ind),'LineWidth',1.5,'LineStyle','-.')
    legend('Experimentally Measured','Smoothing Spline','Power-law fit','$\epsilon_c$','location','Northwest','color','none','Edgecolor','none','Interpreter','Latex','FontName','Arial','Fontsize',22)
    
    Acr(ind) = (mas*ac)/(Elv(ind)*(Epd2(ind)^lam(ind))); % critical area in m^2
    hcr(ind) = ((vi^2)/(2*ac))*((1+lam(ind))/Epd2(ind)); % critical thickness in m

end

%%
Arv_p = [Acr,Are_t];
Arv = 10^(round(log10(min(Arv_p))-1)):10^(round(log10(min(Arv_p))-1)):10^(round(log10(max(Arv_p)))+1);
[mar,iar] = min(abs(Arv-Are_t));

for id1 = 1:length(Arv)

    A_cu = Arv(id1);

    for id2 =1:length(Acr)
        A_rl = Acr(id2);
        
        if A_rl<A_cu
            fac = (A_cu/A_rl).^(1./lam(id2));
        elseif A_rl>A_cu
            fac = A_rl/A_cu;
        elseif A_rl==A_cu
            fac = 1;
        end
        h_rl(id2) = fac.*hcr(id2);
        M_rl(id2) = dens(id2)*h_rl(id2)*A_cu;

    end
    h_rl_mat{id1} = h_rl;
    M_rl_mat{id2} = M_rl;

    [mn,ia] = min(h_rl);
    dat_hmin_abs(id1,1) = ia;
    dat_hmin_abs(id1,2) = mn;

    [mn1,ia1] = min(M_rl);
    dat_Mmin_abs(id1,1) = ia1;
    dat_Mmin_abs(id1,2) = mn1;

end

%% Minimum thickness

figure
yyaxis left
set(gca,'ycolor',[0.64,0.08,0.18])
plot(Arv,dat_hmin_abs(:,2),'Color',[0.64,0.08,0.18],'LineStyle','-.')
hold on
%plot(Acr,hcr,'s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8)
ylabel('$h_{min}\;(m)$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
set(gca,'Xscale','log','LineWidth',1.5)
set(gca,'Yscale','log','LineWidth',1.5)
plot([Arv(iar),Arv(iar)],[1e-3,dat_hmin_abs(iar,2)],'Color','k','LineStyle','-','LineWidth',1.5)
plot(Arv(iar),dat_hmin_abs(iar,2),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
xlim([min(Arv),max(Arv)])
xticks(10.^([log10(min(Arv)):log10(max(Arv))]))
%ylim([1e-3,1])

yyaxis right
set(gca,'ycolor',[0,0,0])
plot(Arv,dat_hmin_abs(:,1),'Color','k','LineStyle','--')
ylabel('Sample Number $(h_{min})$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
set(gca,'Yscale','linear','LineWidth',1.5)
ylim([0,max(dat_hmin_abs(:,1))+1])
xlabel('$A\;(m^2)$','Interpreter','Latex','Rotation',0,'FontName','Arial','color','k','FontSize',24)
ylim([0,1+sz(1)])
yticks([1:sz(1)])
disp(['Minimum thickness is ',num2str(dat_hmin_abs(iar,2)*100),' cm',', achieved using foam ',num2str(dat_hmin_abs(iar,1))])

%% Minimum mass

figure
yyaxis left
set(gca,'ycolor',[0.64,0.08,0.18])
plot(Arv,dat_Mmin_abs(:,2).'./Arv,'Color',[0.64,0.08,0.18],'LineStyle','-.')
hold on
%plot(Acr,hcr.*dens,'s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8)
ylabel('$M_{min}/A\;\;(kg/m^2)$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
set(gca,'Xscale','log','LineWidth',1.5)
set(gca,'Yscale','log','LineWidth',1.5)
plot([Arv(iar),Arv(iar)],[1e-2,dat_Mmin_abs(iar,2)/Arv(iar)],'Color','k','LineStyle','-','LineWidth',1.5)
plot(Arv(iar),dat_Mmin_abs(iar,2)/Arv(iar),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
xlim([min(Arv),max(Arv)])
xticks(10.^([log10(min(Arv)):log10(max(Arv))]))

yyaxis right
set(gca,'ycolor',[0,0,0])
plot(Arv,dat_Mmin_abs(:,1),'Color','k','LineStyle','--')
ylabel('Sample Number $(M_{min})$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
set(gca,'Yscale','linear','LineWidth',1.5)
ylim([0,max(dat_hmin_abs(:,1))+1])
xlabel('$A\;(m^2)$','Interpreter','Latex','Rotation',0,'FontName','Arial','color','k','FontSize',24)
ylim([0,1+sz(1)])
yticks([1:sz(1)])
disp(['Minimum mass per unit area is ',num2str(dat_Mmin_abs(iar,2)/Arv(iar)),' kg/m^2',', achieved using foam ',num2str(dat_Mmin_abs(iar,1))])

%% Geometric constraint

pec = 0:0.01:100;

for ind=1:length(hcr)
  
    [mn,in] = min(abs(Acr(ind)-Arv));
    
    Amin{ind} = Arv(1:in);
    h_Amin{ind} = hcr(ind)*(Acr(ind)./Amin{ind});

    Amax{ind} = Arv(in+1:end);
    h_Amax{ind} = hcr(ind)*((Amax{ind}/Acr(ind)).^(1/lam(ind)));

    %Amin(ind,:) = Acr(ind)*(1./(1+pec));
    %Amax(ind,:) = Acr(ind)*((1+pec).^lam(ind));
    %hrn(ind,:) = hcr(ind)*(1+pec);

    %Avec(ind,:) = [fliplr(Amin(ind,:)),Amax(ind,:)];
    %hvec(ind,:) = [fliplr(hrn(ind,:)),hrn(ind,:)];

    Avec(ind,:) = Arv;
    hvec(ind,:) = [h_Amin{ind},h_Amax{ind}];
    Mvec(ind,:) = hvec(ind,:)*dens(ind);

figure
patch([Avec(ind,1),Avec(ind,:),Avec(ind,end)], [100,hvec(ind,:),100], Colors(ind,:),'FaceAlpha',0.35,'LineStyle','none')
hold on
plot(Amin{ind},h_Amin{ind},'LineStyle','-.','Color','k')
plot(Amax{ind},h_Amax{ind},'LineStyle','-','Color','k')
plot(Acr(ind),hcr(ind),'o','MarkerSize',18,'Color','k','MarkerFaceColor', Colors(ind,:))
set(gca,'Xscale','log','LineWidth',1.5)
set(gca,'Yscale','log','LineWidth',1.5)
xlim([min(Arv),max(Arv)])
xticks(10.^([log10(min(Arv)):log10(max(Arv))]))
ylabel('$h\;(m)$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
xlabel('$A\;(m^2)$','Interpreter','Latex','Rotation',0,'FontName','Arial','color','k','FontSize',24)
ylim([1e-3,1])
box on
end

%%
dife = diff(dat_hmin_abs(:,1));
locsw = [0;find(abs(dife)>0);length(Arv)-1];

figure
plot(Arv,dat_hmin_abs(:,2),'Color',[0.64,0.08,0.18],'LineStyle','-.')
hold on
set(gca,'Xscale','log','LineWidth',1.5)
set(gca,'Yscale','log','LineWidth',1.5)

ULh = 10^(floor(log10(max(dat_hmin_abs(:,2))))+1);
LLh = 10^(floor(log10(min(dat_hmin_abs(:,2))))-2);

for ind=2:length(locsw)
sopt = dat_hmin_abs(locsw(ind),1);  % sample number
lf = locsw(ind-1)+1;
rf = locsw(ind)+1;
patch([Avec(sopt,lf),Avec(sopt,lf:rf),Avec(sopt,rf)], [ULh,hvec(sopt,lf:rf),ULh], Colors(ind,:),'FaceAlpha',0.35,'LineStyle','none')
end

%plot(Acr,hcr,'s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8)
ylabel('$h_{min}\;(m)$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
plot([Arv(iar),Arv(iar)],[LLh,dat_hmin_abs(iar,2)],'Color','k','LineStyle','-','LineWidth',1.5)
plot(Arv(iar),dat_hmin_abs(iar,2),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
xlim([min(Arv),max(Arv)])
xticks(10.^([log10(min(Arv)):log10(max(Arv))]))
%ylim([1e-3,1])
xlabel('$A\;(m^2)$','Interpreter','Latex','Rotation',0,'FontName','Arial','color','k','FontSize',24)
axis tight

%%
dife = diff(dat_Mmin_abs(:,1));
locsw = [0;find(abs(dife)>0);length(Arv)-1];

figure
plot(Arv,dat_Mmin_abs(:,2).'./Arv,'Color',[0.64,0.08,0.18],'LineStyle','-.')
hold on
set(gca,'Xscale','log','LineWidth',1.5)
set(gca,'Yscale','log','LineWidth',1.5)

ULM = 10^(floor(log10(max(dat_Mmin_abs(:,2).'./Arv)))+1);
LLM = 10^(floor(log10(min(dat_Mmin_abs(:,2).'./Arv)))-2);

for ind=2:length(locsw)
sopt = dat_Mmin_abs(locsw(ind),1);  % sample number
lf = locsw(ind-1)+1;
rf = locsw(ind)+1;
patch([Avec(sopt,lf),Avec(sopt,lf:rf),Avec(sopt,rf)], [ULM,Mvec(sopt,lf:rf),ULM], Colors(ind,:),'FaceAlpha',0.35,'LineStyle','none')
end

%plot(Acr,hcr.*dens,'s','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8)
ylabel('$M_{min}/A\;\;(kg/m^2)$','Interpreter','Latex','Rotation',90,'FontName','Arial','FontSize',24)
plot([Arv(iar),Arv(iar)],[LLM,dat_Mmin_abs(iar,2)/Arv(iar)],'Color','k','LineStyle','-','LineWidth',1.5)
plot(Arv(iar),dat_Mmin_abs(iar,2)/Arv(iar),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
xlim([min(Arv),max(Arv)])
xticks(10.^([log10(min(Arv)):log10(max(Arv))]))
xlabel('$A\;(m^2)$','Interpreter','Latex','Rotation',0,'FontName','Arial','color','k','FontSize',24)
axis tight

%%
