function f_n = Gram_Schmidt(x,f,ax,bx)
   
 M = length(f(ax));
 f_n = zeros(length(x),M);
 fx = f(x);
 
 int = 0.001;
 x_seq = linspace(ax,bx,(bx-ax)/int);
 
 f_int = f(x_seq);
 v=f_int(:,1);
 normv=sqrt((v')*v*int);
 f_n(:,1)=fx(:,1)/normv;
 f_n_int = zeros(length(x_seq),M);
 f_n_int(:,1) = f_int(:,1)/normv;
 
 if M>1
 for j=2:M
     proj = (f_int(:,j)')*f_n_int(:,1:j-1).*int;
     proj = repmat(proj,length(x),1);
     projfx = proj.*fx(:,1:j-1);
     %projfx = proj.*f_n(:,1:j-1);
     f_n(:,j)=fx(:,j)-sum(projfx,2);
     
     proj = (f_int(:,j)')*f_n_int(:,1:j-1).*int;
     proj = repmat(proj,length(x_seq),1);
     projfint = proj.*f_int(:,1:j-1);
     %projfint = proj.*f_n_int(:,1:j-1);
     f_n_int(:,j)=f_int(:,j)-sum(projfint,2);
     v= f_n_int(:,j);
     normv = sqrt((v')*v*int);
     
     f_n(:,j)=f_n(:,j)/normv;
     f_n_int(:,j)=f_n_int(:,j)/normv;
     
 end
 end




