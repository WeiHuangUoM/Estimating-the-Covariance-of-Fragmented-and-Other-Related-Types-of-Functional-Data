%tic
%---------------------------------------------------------------
%Smooth the empirical covariance function
%---------------------------------------------------------------

%set spline order

order=4;

extKxl = repmat(ax,1,order);
extKxu = repmat(bx,1,order);
extKyl = repmat(ay,1,order);
extKyu = repmat(by,1,order);


%Set leave one curve out cross-validation criterion

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:SingularMatrix');


if ~exist('K_hat','var')


KK = repmat(Kp,1,length(lambdap));
lambdapp = repmat(lambdap,length(Kp),1);lambdapp = lambdapp(:)';
pp = vertcat(KK,lambdapp);

Cal = arrayfun(@(x)F(pp(1,x),pp(2,x),xs,ys,extKxl,extKxu,extKyl,extKyu,n,order,count,Y,ind,idx1,idx2,t),1:size(pp,2));
paridx = find(Cal==min(Cal));
K = KK(paridx);
lambda1 = lambdapp(paridx);

[N,allknots1,allknots2]= Nform(K,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order);
D=Dform(allknots1,allknots2,n,xs,ys,K,order);
[~,K_hat] = betahat(lambda1,N,D,nw,zs,K,order);


end

%-------------------------------------------------------------
%Fit the smoothed empirical covariance using Fourier period=I
%-------------------------------------------------------------
fsize = floor((n^(1/2)-2)/2);
KK = fsize-2:fsize;
[C,fi] = arrayfun(@(x)LS(KK(x),xs,ys,ax,bx,ay,by,K_hat,int,diag_ind,nw_diag),1:length(KK),'UniformOutput',false);
fi = cell2mat(fi);
ind = find(fi==min(fi));
K_slect = KK(ind);
p = K_slect*2+2;
C = C{ind};
Cv = reshape(C,p^2,1);


if length(xx)>10000
    num = ceil(length(xx)/10000);
    fhat_vec = zeros(1,length(xx));
    for i = 1:num-1
        fhat_vec((i-1)*10000+1:i*10000) = f_hat(xx((i-1)*10000+1:i*10000),yy((i-1)*10000+1:i*10000),ax,bx,ay,by,K_slect,Cv);
    end
    fhat_vec((num-1)*10000+1:length(xx))=f_hat(xx((num-1)*10000+1:length(xx)),yy((num-1)*10000+1:length(xx)),ax,bx,ay,by,K_slect,Cv);
else
    fhat_vec= f_hat(xx,yy,ax,bx,ay,by,K_slect,Cv);
end

fhat_vec = reshape(fhat_vec,1,length(xx));
%-------------------------------------------------------------------------------------------------
%Calculate ISE
%-------------------------------------------------------------------------------------------------
rng(1119*sim);
%for ISE
  J=100;
  n=sample_size;
  delta = fragment_length;
  deltap = delta-0.1;
  l = rand(n,1)*0.2+deltap;
  middle = rand(n,1)*(1-deltap)+deltap/2;
  start = middle-l/2;
  start(start<0)=0;
  endp = middle + l/2;
  endp(endp>1)=1;
  idx1 = ceil(start*(J-1)+1);
  idx2 = floor(endp*(J-1)+1);
  counto = zeros(J);
  for i = 1:n
  
    counto(idx1(i):idx2(i),idx1(i):idx2(i))=counto(idx1(i):idx2(i),idx1(i):idx2(i))+1;
  end

 f_f=reshape(fhat_vec,n1,n2);
 f_f = (f_f+f_f')/2;          
 diff = (f_f - Tru).^2;
 ise_tensor=sum(diff(:)*interval);
 diff_inter = diff;
 diff_inter(counto==0)=nan;
 ise_inter = nansum(diff_inter(:))*interval;           
            
 diff_extend = diff;
 diff_extend(counto~=0)=nan;
 ise_extend = nansum(diff_extend(:))*interval;
            
            

%Set Spline basis matrix
function [N,allknots1,allknots2]=Nform(K,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order)
  p = linspace(0,1, K+2);
  intKx =  quantile(xs,p(2:K+1));
  intKy =  quantile(ys,p(2:K+1));
  
  allknots1=horzcat(extKxl,intKx,extKxu);
  allknots2=horzcat(extKyl,intKy,extKyu);
  
  k=K+order;
  Bdd1 = fnval(spmak(allknots1,eye(k)),xs)';
  Bdd2 = fnval(spmak(allknots2,eye(k)),ys)';
  
  N=zeros(n,k*k);
  
  for i = 1:n
  
    N_temp=kron(Bdd2(i,:),Bdd1(i,:));
    N(i,:)=N_temp;
  end
  
  
end

%Set differentiated basis matrix
function D=Dform(allknots1,allknots2,n,xs,ys,K,order)
  k=K+order;
  Bdd1t = fnval(spmak(allknots1,eye(k)),xs)';
  Bdd2t = fnval(spmak(allknots2,eye(k)),ys)';
  Bdd1q = fnval(fnder(spmak(allknots1,eye(k)),2),xs)';
  Bdd2q = fnval(fnder(spmak(allknots2,eye(k)),2),ys)';
  
  D= zeros(n,size(Bdd1q,2)*size(Bdd2t,2));
  for i = 1:n
    D_temp = kron(Bdd2t(i,:),Bdd1q(i,:))+kron(Bdd2q(i,:),Bdd1t(i,:));
    D(i,:) = D_temp;
  end
end




%coefficients estimation for the tensor product penalized spline
function [coef,estimator] = betahat(lambda1,N,D,nw,zs,K,order)

  lc = size(N,2);
  coef = zeros(1,lc);
  
  col.match=arrayfun(@(x)isequal(N(:,x),zeros(1,size(N,1))),1:lc);
  
  N=N(:,~col.match);
  D=D(:,~col.match);
  
  WB = bsxfun(@times, N, nw);
  BTB = N'*WB;
  
  Dl = D'*D;

 Inverse = pinv(BTB+lambda1*Dl);
 

  for i=1:K+order
    for j=i:K+order
     Inverse((j-1)*(K+order)+i,:)=Inverse((i-1)*(K+order)+j,:);
    end
  end
  
  %Hat=WB*Inverse*WB';
  Cv=Inverse*(WB')*zs;
  estimator=N*Cv;
  coef(~col.match)=Cv;
 
end

%set leave one curve out cv
function cv=F(K,lambda1,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order,count,Y,ind,idx1,idx2,t)

  
  [N,allknots1,allknots2]=Nform(K,xs,ys,extKxl,extKxu,extKyl,extKyu,n,order);
  D=Dform(allknots1,allknots2,n,xs,ys,K,order);
 
  cv=0;
 %for j = 1:length(Yidx)
    %i = Yidx(j);
  for i=1:size(Y,1)
    count_tmp = count;
    count_tmp(idx1(i):idx2(i),idx1(i):idx2(i))=count(idx1(i):idx2(i),idx1(i):idx2(i))-1;
    
    J=length(t);
    Z_tmp = zeros(J,J);
    
     for l = 1:size(Y,1)
      if l==i
          Z_P_temp = zeros(J,J);
      else
         ind_temp=find((idx1>=idx1(l) & idx1<=idx2(l) )| (idx2>=idx1(l) & idx2<=idx2(l)));
         ind_temp(ind_temp==i)=[];
         Y_temp = Y(ind_temp,:);
         Ybarl = nanmean(Y_temp,1);
         Y_tmp = Y(l,:)-Ybarl;
         Z_P_temp=Y_tmp'*Y_tmp;
      end    
        Z_tmp(idx1(l):idx2(l),idx1(l):idx2(l))=Z_tmp(idx1(l):idx2(l),idx1(l):idx2(l))+Z_P_temp(idx1(l):idx2(l),idx1(l):idx2(l));
     end
 q = min(10,ceil(max(count_tmp(:)/20)));
 count_tmp(count_tmp<=q)=0;
 Z_tmp=Z_tmp./count_tmp;  
    
    
    ztmp=reshape(Z_tmp,size(Z_tmp,1)*size(Z_tmp,2),1);
    indt = isnan(ztmp(~ind));
    ztmp=ztmp(~ind);
    ztmp = ztmp(~indt);
    
    Ntmp = N(~indt,:);
    Dtmp = D(~indt,:);
    
    nwtmp = reshape(count_tmp,size(count_tmp,1)*size(count_tmp,2),1);
    nwtmp = nwtmp(~ind);
    nwtmp = nwtmp(~indt);
    
    
    [coef,~]=betahat(lambda1,Ntmp,Dtmp,nwtmp,ztmp,K,order);
    C = reshape(coef,(K+order),(K+order));
    B  = fnval(spmak(allknots1,eye(K+order)),t);
    est_oco = B'*C*B;
 
    Y_tmp = Y(i,:);
    Y_temp = Y((idx1>=idx1(i) & idx1<=idx2(i) )| (idx2>=idx1(i) & idx2<=idx2(i)),:);
    Ybari = nanmean(Y_temp,1);
    Y_temp = Y_tmp-Ybari;
    
    Z_emp = Y_temp'*Y_temp;
    Z_emp = Z_emp(idx1(i):idx2(i),idx1(i):idx2(i));
 
    %cross-validation
    er = (Z_emp-est_oco(idx1(i):idx2(i),idx1(i):idx2(i))).^2;
    error = nansum(er(:));
    cv=cv+error;
  end
  

end




function N=Bform(K,xs,ys,ax,bx,ay,by)

  n=length(xs);
  Bdd1 = Gram_Schmidt(xs,@(x)fourierx(x,K,(bx-ax)),ax,bx);
  Bdd2 = Gram_Schmidt(ys,@(x)fourierx(x,K,(by-ay)),ay,by);
    

  
  N=zeros(n,size(Bdd1,2)*size(Bdd2,2));
  
  for i = 1:n
  
    N_temp=kron(Bdd2(i,:),Bdd1(i,:));
    N(i,:) = N_temp;
  end
end


%unconstraint coefficients estimation for the tensor product Fourier
function coef=fit(WB,BTB,z,p,int)

  
  Inverse = inv(BTB);
  
  
  for i=1:p
    for j=i:p
     Inverse((j-1)*p+i,:)=Inverse((i-1)*p+j,:);
    end
  end
  
 
  coef=Inverse*(WB'*z*int^2);

 
end




%set Sp function  
function [C,fval]=LS(K,xs,ys,ax,bx,ay,by,K_hat,int,diag_ind,nw_diag)

  p=2*K+2;
  
  N=Bform(K,xs,ys,ax,bx,ay,by);
  WB = bsxfun(@times, N, nw_diag);
  BTB = N'*WB*int^2;
  inv(BTB);
  [~, msgid] = lastwarn;
  
  if strcmp(msgid,'MATLAB:singularMatrix')
      C=zeros(p,p);
      fval = Inf;
  elseif strcmp(msgid,'MATLAB:nearlySingularMatrix')
      lastwarn('');
      alpha = eps;
      BTB=BTB+alpha*eye(size(BTB,1));
      inv(BTB);
      [~, msgid] = lastwarn;
      if strcmp(msgid,'MATLAB:nearlySingularMatrix')
         while strcmp(msgid,'MATLAB:nearlySingularMatrix')
         lastwarn('');
         alpha = alpha*100;
         BTB=BTB+alpha*eye(size(BTB,1));
         inv(BTB);
         [~, msgid] = lastwarn;
         end
      end
  
  Cv=fit(WB,BTB,K_hat,p,int);
  
  C0 = reshape(Cv,p,p);
  
  %sprintf('fit')
  
  fun = @(x)Obj(x,p,N,K_hat,int,diag_ind);
 
 
  [Ei_vec,Ei_val] = eig(C0);
  Ei_val(Ei_val<0)=min(Ei_val(Ei_val>0));
  par0 = sqrt(Ei_val)*(Ei_vec)';
  par0 = reshape(par0,p^2,1);
  options= optimoptions(@fminunc,'Algorithm','quasi-Newton','MaxFunctionEvaluations',1000000,'MaxIterations',10000,'SpecifyObjectiveGradient',true,'Display','off');
 
  [par,fval] = fminunc(fun,par0,options);
  V = reshape(par,p,p);
  C = V'*V;

  else
      Cv=fit(WB,BTB,K_hat,p,int);

      C0 = reshape(Cv,p,p);
  
 % sprintf('fit')
  
  fun = @(x)Obj(x,p,N,K_hat,int,diag_ind);
 
 
  [Ei_vec,Ei_val] = eig(C0);
  Ei_val(Ei_val<0)=min(Ei_val(Ei_val>0));
  par0 = sqrt(Ei_val)*(Ei_vec)';
  par0 = reshape(par0,p^2,1);
  options= optimoptions(@fminunc,'Algorithm','quasi-Newton','MaxFunctionEvaluations',1000000,'MaxIterations',10000,'SpecifyObjectiveGradient',true,'Display','off');
 
           [par,fval] = fminunc(fun,par0,options);
            V = reshape(par,p,p); 
            C=V'*V;
  end
end


function [f,g]=Obj(par,p,N,K_hat,int,diag_ind)
  
    V = reshape(par,p,p);
    C = V'*V;
    coef = reshape(C,p^2,1);
    
    int_S = (K_hat'*K_hat - 2*coef'*(N'*K_hat)+ (coef'*N')*N*coef)*int^2;
    int_D = (K_hat(diag_ind)'*K_hat(diag_ind) - 2*coef'*(N(diag_ind,:)'*K_hat(diag_ind))+ (coef'*N(diag_ind,:)')*N(diag_ind,:)*coef)*int; 
    M_D = (K_hat(diag_ind)'*K_hat(diag_ind) - 2*coef'*(N(diag_ind,:)'*K_hat(diag_ind))+ (coef'*N(diag_ind,:)')*N(diag_ind,:)*coef)*int^2; 
    
    f = int_S-M_D+int_D;
    
   if nargout > 1 % supply gradient
    gfc = (-2*(N'*K_hat)+2*(N'*N)*coef)*int^2 -(-2*(N(diag_ind,:)'*K_hat(diag_ind))+2*(N(diag_ind,:)'*N(diag_ind,:))*coef)*int^2+(-2*(N(diag_ind,:)'*K_hat(diag_ind))+2*(N(diag_ind,:)'*N(diag_ind,:))*coef)*int;
    gcx = zeros(p^2,p^2);
    for i = 1:p
        for j=1:p
            
              for s =1:p  
                    gcx((i-1)*p+s,(j-1)*p+i)=par((j-1)*p+s);
                    gcx((j-1)*p+s,(j-1)*p+i)=par((i-1)*p+s);
                    gcx((i-1)*p+s,(i-1)*p+i)=2*par((i-1)*p+s);
                   
                
             end
        end   
    end
    g = gcx*gfc;
    
    end
end

%Constraint estimator
function fhat=f_hat(x,y,ax,bx,ay,by,K,Cv)

  Bs1 = Gram_Schmidt(x,@(x)fourierx(x,K,(bx-ax)),ax,bx);
  Bs2 = Gram_Schmidt(y,@(x)fourierx(x,K,(by-ay)),ay,by);
  
  N=zeros(length(x),size(Bs1,2)*size(Bs2,2));
  
  for i = 1:length(x)
  
    N_temp=kron(Bs2(i,:),Bs1(i,:));
    N(i,:) = N_temp;
  end
  
  fhat = N*Cv;
end