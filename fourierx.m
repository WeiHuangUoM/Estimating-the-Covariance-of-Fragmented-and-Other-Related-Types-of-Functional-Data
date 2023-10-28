%Ferior basis (add x)
function f = fourierx(x,M,period)
 scale = sqrt(2/period);
 f = zeros(length(x),2*M+2);
 f(:,1)=1;
 
 for i=1:M
     f(:,2*i)=cos(2*i*pi*x/period)*scale;
     f(:,2*i+1) = sin(2*i*pi*x/period)*scale;
 end

 f(:,2*M+2)=x;
end
