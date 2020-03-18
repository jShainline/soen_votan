%Shapirosteps
phi=2.067833831e-15;
freq=16e9;
T = 4.7;
pows=[-130,[-34:2:-20,-19:1:-8]];
path = 'E:\DATA\cooldown_march_2019\wafer180327B_Chip24\2019-05-17_ajs\';
cd(path)
%I={};
%V={};
for k=1:length(pows)
    setpower(pows(k))
    pause(20)
    [I_values,V_values]=IVfunction(yokoobj,nanovoltobj);
    I{k}=I_values;
    V{k}=V_values;
    
    save('Shapiro_IV_4_7k_16GHz_ThruMD50.mat','I','V','pows','freq','phi','T')
    
end
setpower(-130)