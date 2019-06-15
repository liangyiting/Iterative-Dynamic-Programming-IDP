function []=erere()
global cp ck ce  cepk gamma etm H betamax betamin;
cp=[25.7 36.6 8.3 2.9 2.1 41.7];%有效 降雨量
ck=[0 0 0 0 0 0];%地下水 不急
ce=[0 0 0 0 0 0];%渗漏梁？？？？
H=[0.5 0.7 0.8 1.0 1.2 1.2];%计划湿润层深度
gamma=[0.0994 0.0406 0.0370 0.2896 0.2809 0.001];
cepk=(ce-cp-ck);
etm=[1 1 1 1 1 1];
betamax=0.9;%突然含水率上线???
betamin=0.1;%土壤含水率下限????
global emax emin;
emax=400;
emin=10;
Q=1200;
dopt=Q/6*ones(6,1);
wopt=dopt;
q0=Q;s0=0;
[qopt,sopt]=fun(1,q0,s0,dopt,wopt);
kssmax=5;%迭代5轮
global N;
N=5;
cq=zeros(N+1,N+1);cs=cq;cf=cq;
cfk=[];
for kss=1:kssmax;
    for i=6:-1:1;
        [cq1,cs1,cf1,cdopt,cwopt]=trans(i,dopt,wopt,qopt,sopt,cq,cs,cf,q0,s0);%由第i+1步骤的状态代价，计算第i步的状态代价(上一次的最优状态轨迹q,s附近）
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
qq=zeros(6-i+1,1);%水存量序列
ss=zeros(6-i+1,1);%含水量序列
q=q0;s=s0;%初始的存水量和含水量
qq(1)=q;ss(1)=s;
cons=[];
for k=i:6;
    dk=dopt(k-i+1);%灌水量
    wk=wopt(k-i+1);%蒸腾
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
%f=@(q,s)interp2(cq,cs,cfq,q,s,'nearest');%二维插值
%q(i),s(i)是在第i个阶段操作之前的状态
di=dopt(i);wi=wopt(i);
dd=50;dw=50;
dq=dd;ds=dd+dw;
if i>1;
    %当前状态的中心是最优轨迹
    q0=qopt(i-1);s0=sopt(i-1);%-cepk(i-1);
end
%dmax=min(q0,di+dd);dmin=max(0,di-dd);
%qmax=q0-dmin;qmin=q0-dmax;
qmax=q0+dq;qmin=max(0,q0-dq);
%wmax=min(emax,wi+dw);wmin=max(emin,wi-dw);%约束 emin<w<emax
%smax=s0+dmax-wmin;smin=s0+dmin-wmax;
smax=s0+ds;smin=max(0,s0-ds);
%获取当前可能的状态--也就是第i-1阶段的结束后的状态（第i阶段的起始状态）
[cq1,cs1]=meshgrid(linspace(qmin,qmax,N),linspace(smin,smax,N));
cdopt=zeros(6,N,N);
cwopt=cdopt;
cf1=zeros(N,N);
for j0=1:N;
    for k=1:N;
        %对正太qi,si寻找最优的d,w以及对应的f和下一个状态的标志,dd，dw是试探步长,需要输入
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

%对每一个状态寻最优的d,w-
function [dopt,wopt,fopt,jopt,kopt]=optim_jk(i,di,wi,qi,si,cf,cq,cs,dd,dw)
global cp ck ce gamma etm H betamax betamin cepk;
global emax emin;
%穷举法
%dd=0.1;dw=0.1;
jmax=5;kmax=5;
cd=-1*ones(2*jmax+1,2*kmax+1);
cw=cd;
cf1=-Inf*cd;
for j0=1:2*jmax+1;
    %提取当前动作d
    dj0=di-dd+j0*dd/jmax;
    %约束1
    if dj0>qi|dj0<0;
        continue;
    end
    %计算下一时刻的q
    qj=qi-dj0;
    cd(j0,:)=dj0;
    for k0=1:2*kmax+1;
        %提取当前动作w
        wj0=wi-dw+k0*dw/kmax;
        %约束2 判断当前动作的可行性
        if wj0<emin|wj0>emax;
            continue;
        end
        %存储动作
        cw(j0,k0)=wj0;
        %计算下一状态的S
        sj=si-cepk(i)+(dj0-wj0);
        %约束3
        Smin=0;Smax=667*H(i)*(betamax-betamin);
        if sj>Smax;
            continue;
        end
        %寻找距离最近的下一个状态-能不能找到最近的点？
        distance=abs(cq-qj)+abs(cs-sj);
        li=min(distance(:));[j1,k1]=find(distance==li);j1=j1(1);k1=k1(1);
        %估算目标函数值
        cf1(j0,k0)=log((wj0/etm(i))^gamma(i))+cf(j1,k1);
    end
end
%寻找使目标函数最大的操作
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
% %可行性约束-保证经过d,w状态转移后步落在第（i+1)步的状态区间内
% %qi-d (min(cq(:)),max(cq(:)));  si+di+pi+ki-ei-w  min(cs(:)),max(cs(:))
% [xopt,fval]=fmincon(@(x)fobj(x(1),x(2)),A,b,[],[],xmin,xmax,opt);
% dopt1=xopt(1);
% wopt1=xopt(2);
% end