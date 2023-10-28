%% 0.Set model
    if CID == 1
     lambda = vertcat(1,0.8,0.3);
      phi = @(x) phi_1(x);
      
    elseif CID == 2
      lambda = vertcat(0.5^(1-1),0.5^(2-1),0.5^(3-1),0.5^(4-1));
      phi = @(x) phi_2(x);
      
    elseif CID == 3
      lambda = vertcat(0.5^(1-1),0.5^(2-1));
      phi = @(x) phi_3(x);
      
    elseif CID == 4
      p=50;
      lambda =(1:p).^(-2);
      phi = @(x,p) phi_4(x,p);  
      
    elseif strcmp(CID,'Demo')
      lambda = vertcat(0.5^(1-1),0.5^(2-1),0.5^(3-1),0.5^(4-1));
      phi = @(x) phi_Demo(x); 
    end

%% 1.Compute the true covariance function

ax = 0; bx = 1; ay=0; by=1; %set the interval of interest
n1=100; n2=100; %number of grid points where the predictions will be computed
interval1 = (bx-ax)/(n1-1);
interval2 = (by-ay)/(n2-1);
interval=interval1*interval2;
xm = ax:interval1:bx;
yn = ay:interval2:by;
xx = repmat(xm,1,n2);yy=repmat(yn,n1,1); yy=yy(:)';
input = vertcat(xx,yy);
tru=arrayfun(@(x) (sum(phi(input(1,x)).*phi(input(2,x)).*lambda)), 1:size(input,2)); %% truth
Tru = reshape(tru,n1,n2);
  
%% 2.Generate data
%% 2.1 Generate full functional data
 rng(1119*sim);
  
  J= 50; 
  grid=linspace(ax,bx,J); %grids of the interval where we can observe data.
  t = (grid-grid(1))/(grid(length(grid))-grid(1));
  int = t(2)-t(1);
  
  n=sample_size;

    kxi = mvnrnd(zeros(n,size(lambda,1)),diag(lambda));
  
    Phi = arrayfun(@(x) phi(grid(1,x)), 1:size(t,2),'UniformOutput', false);
    Phi = cell2mat(Phi);
  
    X = kxi*Phi; %Full functional data
%% 2.2 Generate fragmented functional data 
  delta = fragment_length;
  deltap = delta-0.1;
  %l = rand(n,1)*0.2+deltap;
  l = deltap + 0.1;
  middle = rand(n,1)*(1-deltap)+deltap/2;
  start = middle-l/2;
  start(start<0)=0;
  endp = middle + l/2;
  endp(endp>1)=1;
  
  idx1 = ceil(start*(J-1)+1); %index of the starting point of each fragment
  idx2 = floor(endp*(J-1)+1); %index of the end point of each fragment
  
%% 2.3 Compute the sample covariance function on observed domain  
  count = zeros(J);
  
  for i = 1:n
  
    count(idx1(i):idx2(i),idx1(i):idx2(i))=count(idx1(i):idx2(i),idx1(i):idx2(i))+1;
  end
  
  %Calculate sample convariance on observed domain for Z_P for Panaretos
  Y = zeros(n,J)*nan;
  for  i = 1:n
    Y(i,idx1(i):idx2(i))=X(i,idx1(i):idx2(i));
  end
   Z_P = zeros(J,J);
 for i = 1:n
    
   ind_temp=find((idx1>=idx1(i) & idx1<=idx2(i) )| (idx2>=idx1(i) & idx2<=idx2(i)));
   Y_temp = Y(ind_temp,:);
   Ybari = nanmean(Y_temp,1);
   Y_tmp = Y(i,:)-Ybari;
   Z_P_temp=Y_tmp'*Y_tmp;
  
   Z_P(idx1(i):idx2(i),idx1(i):idx2(i)) =Z_P(idx1(i):idx2(i),idx1(i):idx2(i))+ Z_P_temp(idx1(i):idx2(i),idx1(i):idx2(i));
 end
Z_P = Z_P./count;
Z_test=Z_P;
Z_P(isnan(Z_P))=0;

% Z_ZC: sample covariance on observed domain for ZC method
Ybar=nanmean(Y,1);
Y_c = bsxfun(@minus,Y,Ybar);
Y_c(isnan(Y_c))=0;
Z_ZC=(Y_c'*Y_c)./count;
Z_ZC(isnan(Z_ZC))=0;
  
% Compute sample covariance Z on observed domain for our method
  q =  min(10, ceil(max(count(:))/20)); 
  for i = ceil(J*delta*0.2)+1:J
      for j = 1:ceil(i-J*delta*0.2)
          if count(i,j)<=q
             count(i,j)=nan;
             count(j,i)=nan;
          end
      end
  end
%Calculate empirical convariance for our method
  Z=Z_P.*count./count;
  
%some pre-prossessing of the data for our method  
  z=reshape(Z,size(Z,1)*size(Z,2),1);
  ind = isnan(z);
  zs=z(~ind);
  
  tm = repmat(t,1,J);
  sm = repmat(t,J,1); sm=sm(:)';
  diag_ind = tm==sm;
  xs = tm(~ind);
  ys = sm(~ind);
  diag_ind = diag_ind(~ind);
  
  %set weight matrix
  nw = reshape(count,size(count,1)*size(count,2),1);
  nw = nw(~ind);
 
 
  %set weight for diag_penalize
  nw_diag = nw;
  nw_diag(diag_ind == 0) = 1;
  nw_diag(diag_ind == 1) = 1/int;
  
  
  n=size(xs,2);
  xms = (xm-xm(1))/(xm(length(xm))-xm(1));
  yns=(yn-yn(1))/(yn(length(yn))-yn(1));
  xxs = repmat(xms,1,n2);yys=repmat(yns,n1,1); yys=yys(:)';
  
%Model 1
function y=phi_1(x)

if x<0.5
    y=vertcat((6*x^2-6*x+1)*sqrt(5),log(x+0.5)*sqrt(2),0)*0.5;
else
    y=vertcat( (6*x^2-6*x+1)*sqrt(5),0,(252*x^5-630*x^4+560*x^3-210*x^2+30*x-1)*sqrt(22))*0.5;
end
end
 
%Model 2 (Ruppert)
function y=phi_2(x)
y=vertcat(1, (2*x-1)*sqrt(3), (6*x^2-6*x+1)*sqrt(5), (20*x^3-30*x^2+12*x-1)*sqrt(7));
end

%Model 3
function y=phi_3(x)

y = vertcat(exp(5*(x-0.5))/(1+exp(5*(x-0.5))), (x-0.5)^2/(0.5^2),0,0);
end

%Model 4
function y=phi_4(x,p)

y = zeros(p,length(x));

for i = 1:(p)/2

    y(2*i-1,:)=sqrt(2)*sin(2*i*pi*x);
    y(2*i,:)=sqrt(2)*cos(2*i*pi*x);
   
end

end

%Model 5 (Kneip)
function y=phi_Demo(x)
if x<=0.5
    y=vertcat(sqrt(2)*sin(2*pi*x), sqrt(2)*cos(2*pi*x), 2*sqrt(2)*sin(4*pi*x), 0);
else
    y=vertcat(sqrt(2)*sin(2*pi*x), sqrt(2)*cos(2*pi*x), 0, 2*sqrt(2)*sin(8*pi*x));
end
end
