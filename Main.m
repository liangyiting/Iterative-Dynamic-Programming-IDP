function []=erere()
global cp ck ce  cepk gamma etm H betamax betamin;
cp=[25.7 36.6 8.3 2.9 2.1 41.7];%��Ч ������
ck=[0 0 0 0 0 0];%����ˮ ����
ce=[0 0 0 0 0 0];%��©����������
H=[0.5 0.7 0.8 1.0 1.2 1.2];%�ƻ�ʪ������
gamma=[0.0994 0.0406 0.0370 0.2896 0.2809 0.001];
cepk=(ce-cp-ck);
etm=[1 1 1 1 1 1];
betamax=0.9;%ͻȻ��ˮ������???
betamin=0.1;%������ˮ������????
global emax emin;
emax=400;
emin=10;
Q=1200;
dopt=Q/6*ones(6,1);
wopt=dopt;
q0=Q;s0=0;
[qopt,sopt]=fun(1,q0,s0,dopt,wopt);
kssmax=5;%����5��
global N;
N=5;
cq=zeros(N+1,N+1);cs=cq;cf=cq;
cfk=[];
for kss=1:kssmax;
    for i=6:-1:1;
        [cq1,cs1,cf1,cdopt,cwopt]=trans(i,dopt,wopt,qopt,sopt,cq,cs,cf,q0,s0);%�ɵ�i+1�����״̬���ۣ������i����״̬����(��һ�ε�����״̬�켣q,s������
        cq=cq1;cs=cs1;cf=cf1;
    end
    li=max(cf1(:));[jopt,kopt]=find(cf1==li);jopt=jopt(1);kopt=kopt(1);
    dopt=cdopt(:,jopt,kopt);
    wopt=cwopt(:,jopt,kopt);
    [qopt,sopt,fk]=fun(1,q0,s0,dopt,wopt);
    cfk=[cfk;fk];
end
cfk
figure(1);plot(cfk)
function [qq,ss,fk,cons]=fun(i,q0,s0,dopt,wopt)
global cp ck ce gamma etm H betamax betamin cepk;
global emax emin;
logf=0;
qq=zeros(6-i+1,1);%ˮ��������
ss=zeros(6-i+1,1);%��ˮ������
q=q0;s=s0;%��ʼ�Ĵ�ˮ���ͺ�ˮ��
qq(1)=q;ss(1)=s;
cons=[];
for k=i:6;
    dk=dopt(k-i+1);%��ˮ��
    wk=wopt(k-i+1);%����
    q=q-dk;qq(k-i+1)=q;
    s=s-cepk(k)+dk-wk;ss(k-i+1)=s;    
    logf=logf+log((wk/etm(k))^gamma(k));
    Smax=667*H(k)*(betamax-betamin);
    conk=[s-Smax,-s,-q,-dk];
    cons=[cons;conk];
end
fk=logf;
%fk=exp(logf);

function [cq1,cs1,cf1,cdopt,cwopt]=trans(i,dopt,wopt,qopt,sopt,cq,cs,cf,q0,s0)
global N;
global cp ck ce gamma etm H betamax betamin cepk;
global emax emin;
%f=@(q,s)interp2(cq,cs,cfq,q,s,'nearest');%��ά��ֵ
%q(i),s(i)���ڵ�i���׶β���֮ǰ��״̬
di=dopt(i);wi=wopt(i);
dd=50;dw=50;
dq=dd;ds=dd+dw;
if i>1;
    %��ǰ״̬�����������Ź켣
    q0=qopt(i-1);s0=sopt(i-1);%-cepk(i-1);
end
%dmax=min(q0,di+dd);dmin=max(0,di-dd);
%qmax=q0-dmin;qmin=q0-dmax;
qmax=q0+dq;qmin=max(0,q0-dq);
%wmax=min(emax,wi+dw);wmin=max(emin,wi-dw);%Լ�� emin<w<emax
%smax=s0+dmax-wmin;smin=s0+dmin-wmax;
smax=s0+ds;smin=max(0,s0-ds);
%��ȡ��ǰ���ܵ�״̬--Ҳ���ǵ�i-1�׶εĽ������״̬����i�׶ε���ʼ״̬��
[cq1,cs1]=meshgrid(linspace(qmin,qmax,N),linspace(smin,smax,N));
cdopt=zeros(6,N,N);
cwopt=cdopt;
cf1=zeros(N,N);
for j0=1:N;
    for k=1:N;
        %����̫qi,siѰ�����ŵ�d,w�Լ���Ӧ��f����һ��״̬�ı�־,dd��dw����̽����,��Ҫ����
        [djk,wjk,fjk,jopt,kopt]=optim_jk(i,di,wi,cq1(j0,k),cs1(j0,k),cf,cq,cs,dd,dw);
        cf1(j0,k)=fjk;
        if i<6;
        cdopt(i:end,j0,k)=[djk;cdopt(i+1:end,jopt,kopt)];
        cwopt(i:end,j0,k)=[wjk;cwopt(i+1:end,jopt,kopt)];
        else
             cdopt(i:end,j0,k)=[djk];
        cwopt(i:end,j0,k)=[wjk];
        end
    end
end

%��ÿһ��״̬Ѱ���ŵ�d,w-
function [dopt,wopt,fopt,jopt,kopt]=optim_jk(i,di,wi,qi,si,cf,cq,cs,dd,dw)
global cp ck ce gamma etm H betamax betamin cepk;
global emax emin;
%��ٷ�
%dd=0.1;dw=0.1;
jmax=5;kmax=5;
cd=-1*ones(2*jmax+1,2*kmax+1);
cw=cd;
cf1=-Inf*cd;
for j0=1:2*jmax+1;
    %��ȡ��ǰ����d
    dj0=di-dd+j0*dd/jmax;
    %Լ��1
    if dj0>qi|dj0<0;
        continue;
    end
    %������һʱ�̵�q
    qj=qi-dj0;
    cd(j0,:)=dj0;
    for k0=1:2*kmax+1;
        %��ȡ��ǰ����w
        wj0=wi-dw+k0*dw/kmax;
        %Լ��2 �жϵ�ǰ�����Ŀ�����
        if wj0<emin|wj0>emax;
            continue;
        end
        %�洢����
        cw(j0,k0)=wj0;
        %������һ״̬��S
        sj=si-cepk(i)+(dj0-wj0);
        %Լ��3
        Smin=0;Smax=667*H(i)*(betamax-betamin);
        if sj>Smax;
            continue;
        end
        %Ѱ�Ҿ����������һ��״̬-�ܲ����ҵ�����ĵ㣿
        distance=abs(cq-qj)+abs(cs-sj);
        li=min(distance(:));[j1,k1]=find(distance==li);j1=j1(1);k1=k1(1);
        %����Ŀ�꺯��ֵ
        cf1(j0,k0)=log((wj0/etm(i))^gamma(i))+cf(j1,k1);
    end
end
%Ѱ��ʹĿ�꺯�����Ĳ���
fopt=max(cf1(:));
[jopt,kopt]=find(cf1==fopt);
jopt=jopt(1);kopt=kopt(1);
dopt=cd(jopt,kopt);
wopt=cd(jopt,kopt);
% if 0
% fobj=@(dj,wj)(wi/etm)^gamma(i)*interp2(cq,cs,cfq,qi-dj,si+di+pi+ki-ei-wi,'nearest');
% emin=
% emax=
% %emin<=w<=emax
% %d-w<=..
% qmax_next=max(cq);qmin_next=min(cq);%qmin_next<=qi-d<=qmax_next
% smax_next=max(cs);smin_next=min(cs);%smin<=si+pi+ki-ei+di-wi<=smax
% xmin=[max(0,qi-qmax_next),max([emin,wi-dw])]';xmax=[qi-qmin_next,min([emax,wi+dw])]';
% opt=optimset('display','iter');
% r=
% H=
% A=[1,-1];b=(betamax-betamin)*667*r*H-(cs(i)+cp(i)+ck(i)-ce(i));
% %������Լ��-��֤����d,w״̬ת�ƺ����ڵڣ�i+1)����״̬������
% %qi-d (min(cq(:)),max(cq(:)));  si+di+pi+ki-ei-w  min(cs(:)),max(cs(:))
% [xopt,fval]=fmincon(@(x)fobj(x(1),x(2)),A,b,[],[],xmin,xmax,opt);
% dopt1=xopt(1);
% wopt1=xopt(2);
% end