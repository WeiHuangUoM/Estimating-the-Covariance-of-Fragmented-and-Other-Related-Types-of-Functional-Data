%Set projection matrix
  
  P =zeros(J,J);
  for  pi = 1: J
  for pj  = max(1,(pi+2-floor(J*delta))):min(J,(pi-2+floor(J*delta)) )
  
    P(pi,pj)=1;
  end
  end

  
  
  Z_P(isnan(Z_P))=0;
  
  
 
  [Cp,fi] = arrayfun(@(x) Kf(x,Z_P,P),p:p,'UniformOutput',false);
  fi=cell2mat(fi);
  C = Cp{1};
  
  fhat_vec=arrayfun(@(i)f_hat(xx(i),yy(i),ax,ay,bx,by,J,C),1:length(xx));
  
  
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
  
  function [C,fi]=Kf(k,Z,P)
  
   
    
    [Ei_vec, Ei_val] = eig(Z);
    Ei_vec = fliplr(Ei_vec);
    Ei_val = fliplr(Ei_val);
    Ei_val = flipud(Ei_val);
    Ei_val(Ei_val<0)=0;
    Ei_val = Ei_val(1:k,1:k);
    par =  sqrt(Ei_val)*(Ei_vec(:,1:k)');
    par = reshape(par,k*size(Z,1),1);

    fun = @(x)Obj_P(x,Z,P,k);
    
    options= optimoptions(@fminunc,'Algorithm','quasi-Newton','MaxFunctionEvaluations',1000000,'MaxIterations',10000,'Display','off');
  
    [V,fi] = fminunc(fun,par,options);
    
    V =reshape(V,k,size(Z,1));
    C = V'*V;
   
    
  end
  
   function SE=Obj_P(par,Z,P,k)
      
      d = size(Z,1);
      V = reshape(par,k,d);
      C = V'*V;
      coef = reshape(C,d^2,1);
      
      z = reshape(Z,d^2,1);
      p = reshape(P,d^2,1);
      
      err = (z-p.*coef).^2;
      SE = sum(err);
      
  
   end
   
 
   function f=f_hat(x,y,ax,ay,bx,by,J,C)
  
   f= C(round(((x-ax)/(bx-ax))*(J-1)+1), round(((y-ay)/(by-ay))*(J-1)+1));
   end
  
  