function V=IVthermal(I,Ic,R,T)
%clear all
%xv=linspace(0,2,101);
%T=10;
%Ic=60e-6;
%R=1.6;
%Ic and I in uA;
xv=abs(I)/Ic;
hbar=6.63e-34/(2*pi);
kB=1.38e-23;
e=1.6e-19;
Vc=Ic*R;
gamma= hbar*Ic*1e-6/(e*kB *T);

%%
%fun=@(phi,x,g) exp(-g/2*x*phi).*besseli(0,g*sin(phi/2));
fun=@(phi,x,g)  exp(-g*(x*phi/2-abs(sin(phi/2)))).*besseli(0,g*sin(phi/2),1);
T1=[];
for k=1:length(xv)
    x=xv(k);
    %fun2;
    aux=integral(@(phi) fun(phi,x,gamma),0,2*pi);
    T1=[T1;aux];
end
%%
V=2*Vc/gamma*(1-exp(-pi*gamma*xv)).*(1./T1).*sign(I);

%figure(1)
%plot(V,xv)



