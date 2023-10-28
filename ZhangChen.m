dm = max(idx2-idx1);
a = ceil(0.1*dm);
b = ceil(0.7*dm);
lmax = 1+ceil((J-b)/a);

idl1 = ((1:lmax)-1)*a+1;
idl2 = min(((1:lmax)-1)*a+b,J);

rank = min(length(lambda),b-a);
if CID ==7
    rank =2;
end

Astar = zeros(J,rank);
den = zeros(J,rank);

O1 = eye(rank);

 Sigmal = Z_ZC(idl1(1):idl2(1),idl1(1):idl2(1));
 [Ul,Dl] = eig(Sigmal);
 Ul = fliplr(Ul);
 Dl = rot90(Dl,2);
 Dl(Dl<0)=0;
    
 I1 = idl2(1)-idl1(1)+1;
 sigma2 = (I1-rank)^(-1)*sum(diag(Dl(rank+1:end,rank+1:end)));
    
 Ul = Ul(:,1:rank);
 Dl = Dl(1:rank,1:rank);
    
 D = Dl - sigma2*eye(rank);
 D(D<0) = 0;
 A1 = Ul * D.^(1/2);
 
 Astar(1:I1,:)=A1;
 den(1:I1,:)=den(1:I1,:)+1;
 
for i = 2:lmax
    
    Sigmal = Z_ZC(idl1(i):idl2(i),idl1(i):idl2(i));
    [Ul,Dl] = eig(Sigmal);
    Ul = fliplr(Ul);
    Dl = rot90(Dl,2);
    Dl(Dl<0)=0;
    
    Il = idl2(i)-idl1(i)+1;
    sigma2 = (Il-rank)^(-1)*sum(diag(Dl(rank+1:end,rank+1:end)));
    
    Ul = Ul(:,1:rank);
    Dl = Dl(1:rank,1:rank);
    
    D = Dl - sigma2*eye(rank);
    D(D<0) = 0;
    Al = Ul * D.^(1/2);
    
    USV = (Al(1:(b-a),:))'*A1((a+1):b,:)*O1;
    [U,S,V] = svd(USV);
    
    Ol = U*V';
    
    O1=Ol;
    A1=Al;
    
    Astarl = Al*Ol;
   
    Astar(idl1(i):idl2(i),:)=Astar(idl1(i):idl2(i),:)+ Astarl;
    den(idl1(i):idl2(i),:)=den(idl1(i):idl2(i),:)+1;
end

Atilde = Astar./den;

Sigma0 = Atilde*Atilde';

fhat_vec=arrayfun(@(i)f_hat(xx(i),yy(i),ax,ay,bx,by,J,Sigma0),1:length(xx));

%-------------------------------------------------------------------------------------------------
%Calculate ISE
%-------------------------------------------------------------------------------------------------
  rng(1119*sim);

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
 
 
  function f=f_hat(x,y,ax,ay,bx,by,J,C)
  
   f= C(round(((x-ax)/(bx-ax))*(J-1)+1), round(((y-ay)/(by-ay))*(J-1)+1));
   end
  