
A=1;
m=0;
n=0;
x=linspace(0,1.4,100);
y=linspace(0,1.0,100);
coord = [];
f0=130.8;
ImpTrans=[];
f=[];
MN=[];
Zmn=zeros(100);
   for m=0:50
       for n=0:50
            f(m+1,n+1)=0.453*5447.25*0.003*((((m+1)^2)/(1.4^2))+((n+1)^2))
            if f(m+1,n+1)<=f0
                A=0.0004*(f(m+1,n+1))/2*f0;
            end
            if f(m+1,n+1)>f0
                A=0.0004*f0/(2*f(m+1,n+1));
            end
            for k=1:100
                for j=1:100
                    Z(k,j)=A*((sin((m+1)*3.14*x(k)/1.4)*sin((n+1)*3.14*y(j)/1)));
                end
            end
           Zmn=Zmn+Z;
       end
   end

surf(x,y,Zmn);