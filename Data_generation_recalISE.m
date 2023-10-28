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
  

%--------------------------------------------------------------------------
%Generate data
%--------------------------------------------------------------------------
  rng(1119*sim);
  
  J= 100;
  
  t=linspace(ax,bx,J);
  
  int = t(2)-t(1);
  
  n=sample_size;
   kxi = mvnrnd(zeros(n,size(lambda,1)),diag(lambda));
  
    Phi = arrayfun(@(x) phi(t(1,x)), 1:size(t,2),'UniformOutput', false);
    Phi = cell2mat(Phi);
  
    X = kxi*Phi;
    

%Fragmentisation
  delta = fragment_length;
  deltap = delta-0.1;
  %l = rand(n,1)*0.2+deltap;
  l = deltap+0.1;
  middle = rand(n,1)*(1-deltap)+deltap/2;
  start = middle-l/2;
  start(start<0)=0;
  endp = middle + l/2;
  endp(endp>1)=1;
  
  idx1 = ceil(start*(J-1)+1);
  idx2 = floor(endp*(J-1)+1);
  
  
  count = zeros(J);
  
  for i = 1:n
  
    count(idx1(i):idx2(i),idx1(i):idx2(i))=count(idx1(i):idx2(i),idx1(i):idx2(i))+1;
  end

  
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
