%%script for random distrubution of NPs in total LN pop
%% initialize
clear
h=waitbar(0,'Model running');
n_iteration=10000;

%% loop to get statistics
for i=1:n_iteration

    %% randomly determine the number of each cells based on mean and sd data
    n_tot=round(347 + (randn*9.53));
    
    n_mip=round(146.4167 + (randn*38.4755));
    if n_mip>n_tot, n_mip=n_tot;end;
    if n_mip<0, n_mip=0;end;
    
    n_atr=round(92.03125 + (randn*19.74958));
    if n_atr>n_tot, n_atr=n_tot;end;
    if n_atr<0, n_atr=0;end;
    
    n_tkk=round(9.166 + (randn*2.57));
    if n_tkk>n_tot, n_tkk=n_tot;end;
    if n_tkk<0, n_tkk=0;end;
    
    n_fmrf=round(41.3889 + (randn*7.904));
    if n_fmrf>n_tot, n_fmrf=n_tot;end;
    if n_fmrf<0, n_fmrf=0;end;

    %% randomly assign the different cells "positiveness" to cell index
    bigM=zeros(n_tot,7);
    p=randperm(n_tot);
    bigM(p(1:n_mip),3)=1;
    p=randperm(n_tot);
    bigM(p(1:n_atr),4)=1;
    p=randperm(n_tot);
    bigM(p(1:n_tkk),5)=1;
    p=randperm(n_tot);
    bigM(p(1:n_fmrf),6)=1;

    %% create a "neuropeptidergic" column
    bigM(:,7)=sum(bigM(:,2:6),2);

    ind=find(bigM(:,7)>1);
    bigM(ind,7)=1;

    %% count number of cells that express pairs of NT|NP
    
    n_tkk_n_fmrf(i)=length(find(bigM(:,5)>0 & bigM(:,6)>0));%%TKK 
    n_tkk_n_atr(i)=length(find(bigM(:,5)>0 & bigM(:,4)>0));
    n_tkk_n_mip(i)=length(find(bigM(:,5)>0 & bigM(:,3)>0));%%TKK

    n_mip_n_atr(i)=length(find(bigM(:,3)>0 & bigM(:,4)>0));%%MIP
    n_mip_n_tkk(i)=length(find(bigM(:,3)>0 & bigM(:,5)>0));
    n_mip_n_fmrf(i)=length(find(bigM(:,3)>0 & bigM(:,6)>0));  
    %%MIP
    
    n_atr_n_tkk(i)=length(find(bigM(:,4)>0 & bigM(:,5)>0));%%ATR
	n_atr_n_fmrf(i)=length(find(bigM(:,4)>0 & bigM(:,6)>0));
    n_atr_n_mip(i)=length(find(bigM(:,4)>0 & bigM(:,3)>0));%%ATR
    
    n_fmrf_n_mip(i)=length(find(bigM(:,6)>0 & bigM(:,3)>0));%%FMRF
    n_fmrf_n_atr(i)=length(find(bigM(:,6)>0 & bigM(:,4)>0));
    n_fmrf_n_tkk(i)=length(find(bigM(:,6)>0 & bigM(:,5)>0));%%FMRF

       
    %% calculate percentage of X cells that express Y
    
    p_mip_thatHas_fmrf(i)=n_mip_n_fmrf(i)/n_mip;%%for MIP
    p_mip_thatHas_atr(i)=n_mip_n_atr(i)/n_mip;
    p_mip_thatHas_tkk(i)=n_mip_n_tkk(i)/n_mip;%%for MIP
    
    p_tkk_thatHas_fmrf(i)=n_tkk_n_fmrf(i)/n_tkk;%%for TKK
    p_tkk_thatHas_atr(i)=n_tkk_n_atr(i)/n_tkk;
    p_tkk_thatHas_mip(i)=n_tkk_n_mip(i)/n_tkk;%%for TKK
    
    p_atr_thatHas_fmrf(i)=n_atr_n_fmrf(i)/n_atr;%%for ATR
    p_atr_thatHas_tkk(i)=n_atr_n_tkk(i)/n_atr;
    p_atr_thatHas_mip(i)=n_atr_n_mip(i)/n_atr;%%for ATR
    
    
    p_fmrf_thatHas_tkk(i)=n_fmrf_n_tkk(i)/n_fmrf;%%for FMRF
    p_fmrf_thatHas_atr(i)=n_fmrf_n_atr(i)/n_fmrf;
    p_fmrf_thatHas_mip(i)=n_fmrf_n_mip(i)/n_fmrf;%%for FMRF
    
   waitbar(i/n_iteration,h);
    
end; 

%% plotting numbers of cells that co-express each pairwise combination of neuropeptides

figure;hist(n_mip_n_atr,[0.5:max(n_mip_n_atr)]);title('n_mip_n_atr');%%MIP
figure;hist(n_mip_n_tkk,[0.5:max(n_mip_n_tkk)]);title('n_mip_n_tkk');
figure;hist(n_mip_n_fmrf,[0.5:max(n_mip_n_fmrf)]);title('n_mip_n_fmrf');

figure;hist(n_tkk_n_fmrf,[0.5:max(n_tkk_n_fmrf)]);title('n_tkk_n_fmrf');%%TKK
figure;hist(n_tkk_n_atr,[0.5:max(n_tkk_n_atr)]);title('n_tkk_n_atr');
figure;hist(n_tkk_n_mip,[0.5:max(n_tkk_n_mip)]);title('n_tkk_n_mip');%%TKK

figure;hist(n_atr_n_tkk,[0.5:max(n_atr_n_tkk)]);title('n_atr_n_tkk');%%ATR
figure;hist(n_atr_n_fmrf,[0.5:max(n_atr_n_fmrf)]);title('n_atr_n_fmrf');
figure;hist(n_atr_n_mip,[0.5:max(n_atr_n_mip)]);title('n_atr_n_mip');%%ATR

figure;hist(n_fmrf_n_tkk,[0.5:max(n_fmrf_n_tkk)]);title('n_fmrf_n_tkk');%%FMRF
figure;hist(n_fmrf_n_atr,[0.5:max(n_fmrf_n_atr)]);title('n_fmrf_n_atr');
figure;hist(n_fmrf_n_mip,[0.5:max(n_fmrf_n_mip)]);title('n_fmrf_n_mip');%%FMRF


%% plotting precentage of coexpressing
%% plotting percentage of coexpressing

figure;hist(p_mip_thatHas_fmrf.*100,[0.5:99.5]);title('% MIP cells expressing FMRF');%%MIP
axis([0 100 0 10000])
figure;hist(p_mip_thatHas_atr.*100,[0.5:99.5]);title('% MIP cells expressing ATR');
axis([0 100 0 10000])
figure;hist(p_mip_thatHas_tkk.*100,[0.5:99.5]);title('% MIP cells expressing TKK');%%MIP
axis([0 100 0 10000])

figure;hist(p_atr_thatHas_fmrf.*100,[0.5:99.5]);title('% ATR cells expressing FMRF');%%ATR
figure;hist(p_atr_thatHas_mip.*100,[0.5:99.5]);title('% ATR cells expressing MIP');
axis([0 100 0 10000])
figure;hist(p_atr_thatHas_tkk.*100,[0.5:99.5]);title('% ATR cells expressing TKK');%%ATR
axis([0 100 0 10000])

figure;hist(p_tkk_thatHas_fmrf.*100,[0.5:99.5]);title('% TKK cells expressing FMRF');%%TKK
axis([0 100 0 10000])
figure;hist(p_tkk_thatHas_mip.*100,[0.5:99.5]);title('% TKK cells expressing MIP');
axis([0 100 0 10000])
figure;hist(p_tkk_thatHas_atr.*100,[0.5:99.5]);title('% TKK cells expressing ATR');%%TKK


figure;hist(p_fmrf_thatHas_tkk.*100,[0.5:99.5]);title('% FMRF cells expressing TKK');%%FMRF
axis([0 100 0 10000])
figure;hist(p_fmrf_thatHas_mip.*100,[0.5:99.5]);title('% FMRF cells expressing MIP');
axis([0 100 0 10000])
figure;hist(p_fmrf_thatHas_atr.*100,[0.5:99.5]);title('% FMRF cells expressing ATR');%%FMRF
axis([0 100 0 10000])

close(h);

save random_probability_script_lizbinski.mat
%% the averaging script can be used to determine the average % co-expression and standard deviation of every pairwise 
%%co-expression relationship