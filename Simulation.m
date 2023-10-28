CID='Demo';

sample_size=50;
fragment_length=0.4;
par_m_Pen5040=zeros(100,6); %Parameters
Est_m_Pen5040=cell(100,1); %covriance
par_P5040=zeros(100,4);
Est_P5040=cell(100,1);
par_ZC5040=zeros(100,3);
Est_ZC5040=cell(100,1);

for sim =1:100
  
Data_generation
p=4;
Panaretos
par_P5040(sim,:) =  [ise_tensor,ise_inter,ise_extend,p];
Est_P5040{sim}= f_f;
sprintf('DP: sim %d, %f, %f', sim, ise_inter,ise_extend) 
 
Data_generation
clear K_hat
Kp = 5:15;
lambdap=exp(-15);
Penalize_diagonal
par_m_Pen5040(sim,:) =  [ise_tensor,ise_inter,ise_extend, K,lambda1,p];
Est_m_Pen5040{sim}= f_f;
sprintf('Fourier: sim %d, %f, %f, %d, %f', sim, ise_inter,ise_extend, K, lambda1)


Data_generation
ZhangChen
par_ZC5040(sim,:) =  [ise_tensor,ise_inter,ise_extend];
Est_ZC5040{sim}= f_f;
sprintf('ZC: sim %d, %f, %f', sim, ise_inter,ise_extend) 
end

save('par_m_Pen5040','par_m_Pen5040')
save('Est_m_Pen5040','Est_m_Pen5040')
save('par_P5040','par_P5040')
save('Est_P5040','Est_P5040')
save('par_ZC5040','par_ZC5040')
save('Est_ZC5040','Est_ZC5040')

