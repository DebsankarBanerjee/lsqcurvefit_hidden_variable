clear all
clc 

%------------------------------------------Start taking data-------------------------------------------
name='m4c2';
global ndstart ndend nd Rmax
ndstart=1;
ndend=5;
nd=[ndstart:ndend]; 

names=[];
Ravg=zeros(numel(nd),1);
nt=zeros(numel(nd),1);

for i=1:numel(nd)
char = strcat('data/',name,'_droplet_',num2str(nd(i)),'.txt');
str = string(char);
names=[names;str];
end

global Tnuc
Tnuc=zeros(numel(nd),1);

for i=1:numel(nd)
dropletname = names(i,:); 
d = dropletname(1); 
dat=readmatrix(dropletname(1));
nt=numel(dat(:,1));
tspan(1:nt,i)=dat(:,1)./60;	%% input time in sec - convert to min
Rdata(1:nt,i)=dat(:,2)./2;	%% input was diameter - convert to radius
Tnuc(i)=dat(1,1)./60;		%% nuc time in min
Ravg(i)=mean(Rdata(nt-5:nt,i));

end

tdata=tspan(:,1);		%% one array of time points - same for all droplets as they are measured from the same frame

%------------------------------ rearrange to put zeros in the beginning ---------------------
for i=2:numel(nd)
tdum=tspan(:,i);
tdum(tdum==0)=[];
ndf=numel(tspan(:,1))-numel(tdum);
adum=zeros(ndf,1);
tdum=[adum; tdum];
tspan(:,i)=tdum;		%% not using it in the current version

rdum=Rdata(:,i);
rdum(rdum==0)=[];
rdum=[adum; rdum];
Rdata(:,i)=rdum;
end
%---------------------------------------End taking data-------------------------------------

%---------------------------------------Try to fit -----------------------------------------

% ------------ guess parameters (4 + droplet): a, b, alpha, beta(i) and R0 ----------

Rmax = 0.15;
beta_guess = 0.001* ones(1,numel(nd)); 
beta_stiff_guess = 0.01;
parguess_init = [0.015 0.006 0.01 50]; 						% nonlinear fits need initial guesses
parguess = horzcat(parguess_init,beta_guess,beta_stiff_guess);
lb = zeros(1, 4+ numel(nd)+1);
ub_init = [100 100 100 10000];		% Inf
ub_beta = 1.0*ones(1,numel(nd));
ub_beta_st = 2;
 					 
ub = horzcat (ub_init, ub_beta,ub_beta_st); 
options = optimoptions('lsqcurvefit', 'PlotFcn', 'optimplotresnorm', 'StepTolerance', 1E-9);

%% options: 'Algorithm','levenberg-marquardt', 'MaxIterations',10,
%% more at: https://www.mathworks.com/help/optim/ug/lsqcurvefit.html

[pars, resnorm]  = lsqcurvefit(@myfunc,parguess,tdata,Rdata,lb,ub,options)

display(pars(1:4))
display(pars(4+1:4+numel(nd)))


C = {'k',[0.9290, 0.6940, 0.1250],'r','g',[0.4940, 0.1840, 0.5560],[.5 .6 .7]};

figure(5)
hold on
for i = 1:numel(nd)
    plot(tdata,Rdata(:,i),'o')
end 
set(gca,'FontSize',18)
pbaspect([1 1 1])
title('droplets','Interpreter', 'none');
hold on

tfit = linspace(tspan(1),tspan(end),size(tspan,1));
fit = myfunc(pars,tfit); 
for i = 1:numel(nd)
plot(tfit,fit(:,i),'-','linewidth',2,'color',C{i})
end
hold off

if (resnorm<2)
%---------- write fit results --------
fid2=fopen(['fit/' 'fitparam' '.txt'],'a+');
fprintf(fid2, '%s\n', names(1) );
fprintf(fid2, '%s\n', "-------- guess values -------" );
for i=1:numel(parguess)
fprintf(fid2, '%f\n', parguess(i) );
end
fprintf(fid2, '%s\n', "-------- fit values -------" );
fprintf(fid2, '%s\n', "Fitted using Model 4" );
fprintf(fid2, '%s %f\n', "Resnorm=",resnorm );
for i=1:numel(parguess)
fprintf(fid2, '%f\n', pars(i) );
end
%------------------
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = myfunc(pars,tspan)
        global Tnuc ndstart ndend nd Rmax
        k = pars;
        r0=0.005; 
            function dRdt = ode(t,R)
		        er=1e-4;
                nn=numel(nd);
                % this is the ODE we are fitting to
                dRdt = zeros(nn+1,1);
                for i = 1:nn;
                    dRdt(i) = heaviside(t-Tnuc(i)) * ( k(1)./(R(i)+er) - (k(2)./(R(i)+er)).*(sum(R(1:nn).^3)) - (k(2)*k(4))./(R(i)+er).*(R(end).^3) - k(3).*exp( k(4+i) * (R(i)/(Rmax+R(i))) )./(R(i)+er) );
                end
                dRdt(end) = k(1)/(R(end)+er) - (k(2)/(R(end)+er))*(sum(R(1:nn).^3)) - (k(2)*k(4))*R(end)*R(end) - k(3)*exp( k(4+nn+1) * (R(end)/(Rmax+R(end))) )/(R(end)+er);
		%fprintf(fid3, '%f\n', dRdt );
            end

        R_mat = r0*ones(numel(nd)+1,1);  
        [t,R] = ode15s(@ode,tspan,R_mat);
        dum=R(:,1:end-1);
	R=dum;
        output = R; 
    
end % myfunc





