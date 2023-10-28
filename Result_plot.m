    

    axint = 0;
    bxint = 1;

    ayint = 0;
    byint = 1;

    n1=100;n2=100;
    interval1 = (bxint-axint)/(n1-1);
    interval2 = (byint-ayint)/(n2-1);
    interval=interval1*interval2;
    xm = linspace(axint,bxint,n1);
    yn=linspace(ayint,byint,n2);

    
sP5020 = sort(par_P5020(:,3));
idx = find(par_P5020(:,3)==sP5020(50));
f_f=Est_P5020{idx};
filename = sprintf('SL-5020-Kneip1-Panaretos.pdf');
sample_size=50;
fragment_length=0.2;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sP5050 = sort(par_P5050(:,3));
idx = find(par_P5050(:,3)==sP5050(51));
f_f=Est_P5050{idx};
filename = sprintf('SL-5050-Kneip1-Panaretos.pdf');
sample_size=50;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sP10020 = sort(par_P10020(:,3));
idx = find(par_P10020(:,3)==sP10020(50));
f_f=Est_P10020{idx};
filename = sprintf('SL-10020-Kneip1-Panaretos.pdf');
sample_size=100;
fragment_length=0.2;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sP10050 = sort(par_P10050(:,3));
idx = find(par_P10050(:,3)==sP10050(50));
f_f=Est_P10050{idx};

filename = sprintf('SL-10050-Kneip1-Panaretos.pdf');
sample_size=100;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 



sP50020 = sort(par_P50020(:,3));
idx = find(par_P50020(:,3)==sP50020(50));
f_f=Est_P50020{idx};
filename = sprintf('SL-50020-Kneip1-Panaretos.pdf');
sample_size=500;
fragment_length=0.2;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sP50050 = sort(par_P50050(:,3));
idx = find(par_P50050(:,3)==sP50050(50));
f_f=Est_P50050{idx};
filename = sprintf('SL-50050-Kneip1-Panaretos.pdf');
sample_size=500;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sm5020 = sort(par_m_Pen5020(:,3));
idx = find(par_m_Pen5020(:,3)==sm5020(50));
f_f=Est_m_Pen5020{idx};
filename = sprintf('SL-5020-Kneip1-fx.pdf');
sample_size=50;
fragment_length=0.2;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-2,3])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sm5050 = sort(par_m_Pen5050(:,3));
idx = find(par_m_Pen5050(:,3)==sm5050(51));
f_f=Est_m_Pen5050{idx};
filename = sprintf('SL-5050-Kneip1-fx.pdf');
sample_size=50;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sm10020 = sort(par_m_Pen10020(:,3));
idx = find(par_m_Pen10020(:,3)==sm10020(50));
f_f=Est_m_Pen10020{idx};
filename = sprintf('SL-10020-Kneip1-fx.pdf');
sample_size=100;
fragment_length=0.2;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sm10050 = sort(par_m_Pen10050(:,3));
idx = find(par_m_Pen10050(:,3)==sm10050(50));
f_f=Est_m_Pen10050{idx};
filename = sprintf('SL-10050-Kneip1-fx.pdf');
sample_size=100;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sm50020 = sort(par_m_Pen50020(:,3));
idx = find(par_m_Pen50020(:,3)==sm50020(5));
f_f=Est_m_Pen50020{idx};
filename = sprintf('SL-50020-Kneip1-fx.pdf');
sample_size=500;
fragment_length=0.2;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sm50050 = sort(par_m_Pen50050(:,2));
idx = find(par_m_Pen50050(:,2)==sm50050(50));
f_f=Est_m_Pen50050{idx};
filename = sprintf('SL-50050-Kneip1-fx.pdf');
sample_size=500;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sZC5050 = sort(par_ZC5050(:,3));
idx = find(par_ZC5050(:,3)==sZC5050(51));
f_f=Est_ZC5050{idx};
filename = sprintf('SL-5050-Kneip1-ZC.pdf');
sample_size=50;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 

sZC50050 = sort(par_ZC50050(:,3));
idx = find(par_ZC50050(:,3)==sZC50050(51));
f_f=Est_ZC50050{idx};
filename = sprintf('SL-50050-Kneip1-ZC.pdf');
sample_size=500;
fragment_length=0.5;
sim=idx;
Data_generation_recalISE
f_f_B = f_f;
f_f_B(count==0)=nan;
f_f_off = f_f;
f_f_off(count~=0)=nan;
surf(xm,yn,f_f_B,'FaceAlpha',0.1);
hold on
surf(xm,yn,f_f_off);
hold off
zlim([-5,5])
xlabel('t','FontSize',20)
ylabel('s','FontSize',20)
zlabel('K(s,t)','FontSize',20)
set(gca, 'fontsize', 20); set(gcf, 'PaperPosition', [0 0 20 20]);  set(gcf, 'PaperSize', [20 20]);  print(gcf,'-dpdf',filename, '-opengl') 
