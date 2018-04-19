%% initialize
clear
h=waitbar(0,'Model running');
n_iteration=10000;

%% loop to get statistics
for i=1:n_iteration

    %% randomly determine the number of each cells based on mean and sd data
    n_tot=round(153 + (randn*6.609));
    
    n_mip=round(146.4167 + (randn*38.47));
    if n_mip>n_tot, n_mip=n_tot;end;
    if n_mip<0, n_mip=0;end;
    
    n_atr=round(92.03 + (randn*19.749));
    if n_atr>n_tot, n_atr=n_tot;end;
    if n_atr<0, n_atr=0;end;
    if n_atr>n_mip, n_atr=n_mip;end;
    
    n_tkk=round(9.166 + (randn*2.57));
    if n_tkk>n_tot, n_tkk=n_tot;end;
    if n_tkk<0, n_tkk=0;end;
    
    n_fmrf=round(41.38 + (randn*7.9));
    if n_fmrf>n_tot, n_fmrf=n_tot;end;
    if n_fmrf<0, n_fmrf=0;end;

    %% randomly assign the different cells "positiveness" to cell index
    bigM=zeros(n_tot,7);
    
    p_tkk=randperm(n_tot);
    bigM(p_tkk(1:n_tkk),5)=1;
      
    p_fmrf=randperm(n_tot);
    bigM(p_fmrf(1:n_fmrf),6)=1;
   
    p_mip=randperm(n_tot);
    bigM(p_mip(1:n_mip),3)=1;
    
  %% assumption of 43.59% of FMRF cells are ATRergic (this can be changed to reflect any co-expression relationship 
  %%that you want to specify as a 'rule'
    ind_fmrf=p_fmrf(1:n_fmrf);
    ind_sh=randperm(length(ind_fmrf));
    ind_fmrf_atr=ind_fmrf(ind_sh(1:round(length(ind_sh)*.4359)));
    bigM(ind_fmrf_atr,4)=1;
    
    % distribute the rest of the ATRergic cells among the non_FMRF cells
    indXX=find(bigM(:,4)==1);
    n_atr_rest=n_atr-length(indXX);
    
    ind=find(bigM(:,6)==0);
    p=randperm(length(ind));
    
    try,
        bigM(ind(p(1:n_atr_rest)),4)=1;
    catch
        if n_atr_rest>length(p),n_atr_rest=length(p);end
        if n_atr_rest<length(p),n_atr_rest=0;end
        bigM(ind(p(1:n_atr_rest)),4)=1;
    end
    
    tot_atr_cells(i)=length(find(bigM(:,4)==1));
    
    
    
    %% create a "neuropeptidergic" column
    bigM(:,7)=sum(bigM(:,2:6),2);

    ind=find(bigM(:,7)>1);
    bigM(ind,7)=1;

    %% count number of cells that express pairs of NT|NP
    
    n_tkk_n_fmrf(i)=length(find(bigM(:,5)>0 & bigM(:,6)>0));%%TKK 
    if n_tkk_n_fmrf(i)> n_tkk, n_tkk_n_fmrf(i)=n_tkk; end
    n_tkk_n_atr(i)=length(find(bigM(:,5)>0 & bigM(:,4)>0));
    if n_tkk_n_atr(i)> n_tkk, n_tkk_n_atr(i)=n_tkk; end
    n_tkk_n_mip(i)=length(find(bigM(:,5)>0 & bigM(:,3)>0));%%TKK
    if n_tkk_n_mip(i)> n_tkk, n_tkk_n_mip(i)=n_tkk; end

    n_mip_n_atr(i)=length(find(bigM(:,3)>0 & bigM(:,4)>0));%%MIP
    if n_mip_n_atr(i)> n_mip, n_mip_n_atr(i)=n_mip;end
    n_mip_n_tkk(i)=length(find(bigM(:,3)>0 & bigM(:,5)>0));
    if n_mip_n_tkk(i)> n_mip, n_mip_n_tkk(i)=n_mip;end
    n_mip_n_fmrf(i)=length(find(bigM(:,3)>0 & bigM(:,6)>0)); 
    if n_mip_n_fmrf(i)> n_mip, n_mip_n_fmrf(i)=n_mip;end %%MIP
    
    
    %%ATR
    n_atr_n_tkk(i)=length(find(bigM(:,4)>0 & bigM(:,5)>0));
    if n_atr_n_tkk(i)> n_atr, n_atr_n_tkk(i)=n_atr;end
	n_atr_n_fmrf(i)=length(find(bigM(:,4)>0 & bigM(:,6)>0));
    if n_atr_n_fmrf(i)> n_atr, n_atr_n_fmrf(i)=n_atr;end
    n_atr_n_mip(i)=length(find(bigM(:,4)>0 & bigM(:,3)>0));
    if n_atr_n_mip(i)> n_atr, n_atr_n_mip(i)=n_atr;end
    %%ATR
    
    %%FMRF
    n_fmrf_n_mip(i)=length(find(bigM(:,6)>0 & bigM(:,3)>0));
    if n_fmrf_n_mip(i) > n_fmrf, n_fmrf_n_mip(i)=n_fmrf; end
    n_fmrf_n_atr(i)=length(find(bigM(:,6)>0 & bigM(:,4)>0));
    if n_fmrf_n_atr(i) > n_fmrf, n_fmrf_n_atr(i)=n_fmrf; end
    n_fmrf_n_tkk(i)=length(find(bigM(:,6)>0 & bigM(:,5)>0));%%FMRF
    if n_fmrf_n_tkk(i) > n_fmrf, n_fmrf_n_tkk(i)=n_fmrf; end
    
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

%% plotting n amount of cells that co-express each pairwise combination

%%figure;hist(n_mip_n_atr,[0.5:max(n_mip_n_atr)]);title('n_mip_n_atr');%%MIP
%%figure;hist(n_mip_n_tkk,[0.5:max(n_mip_n_tkk)]);title('n_mip_n_tkk');
%%figure;hist(n_mip_n_fmrf,[0.5:max(n_mip_n_fmrf)]);title('n_mip_n_fmrf');
%%figure;hist(n_mip_n_ast,[0.5:max(n_mip_n_ast)]);title('n_mip_n_ast');%%MIP

%%figure;hist(n_ast_n_atr,[0.5:max(n_ast_n_atr)]);title('n_ast_n_atr');%%AST
%%figure;hist(n_ast_n_tkk,[0.5:max(n_ast_n_tkk)]);title('n_ast_n_tkk');
%%figure;hist(n_ast_n_fmrf,[0.5:max(n_ast_n_fmrf)]);title('n_ast_n_fmrf');
%%figure;hist(n_ast_n_mip,[0.5:max(n_ast_n_mip)]);title('n_ast_n_mip');%%AST

%%figure;hist(n_tkk_n_fmrf,[0.5:max(n_tkk_n_fmrf)]);title('n_tkk_n_fmrf');%%TKK
%%figure;hist(n_tkk_n_atr,[0.5:max(n_tkk_n_atr)]);title('n_tkk_n_atr');
%%figure;hist(n_tkk_n_ast,[0.5:max(n_tkk_n_ast)]);title('n_tkk_n_ast');
%%figure;hist(n_tkk_n_mip,[0.5:max(n_tkk_n_mip)]);title('n_tkk_n_mip');%%TKK

%%figure;hist(n_atr_n_tkk,[0.5:max(n_atr_n_tkk)]);title('n_atr_n_tkk');%%ATR
%%figure;hist(n_atr_n_fmrf,[0.5:max(n_atr_n_fmrf)]);title('n_atr_n_fmrf');
%%figure;hist(n_atr_n_ast,[0.5:max(n_atr_n_ast)]);title('n_atr_n_ast');
%%figure;hist(n_atr_n_mip,[0.5:max(n_atr_n_mip)]);title('n_atr_n_mip');%%ATR

%%figure;hist(n_fmrf_n_tkk,[0.5:max(n_fmrf_n_tkk)]);title('n_fmrf_n_tkk');%%FMRF
%%figure;hist(n_fmrf_n_atr,[0.5:max(n_fmrf_n_atr)]);title('n_fmrf_n_atr');
%%figure;hist(n_fmrf_n_ast,[0.5:max(n_fmrf_n_ast)]);title('n_fmrf_n_ast');
%%figure;hist(n_fmrf_n_mip,[0.5:max(n_fmrf_n_mip)]);title('n_fmrf_n_mip');%%FMRF

%% plotting percentage of cells that co-express each pairwise combination of neuropeptides
figure;hist(p_mip_thatHas_fmrf.*100,[0.5:99.5]);title('% MIP cells expressing FMRF');%%MIP
figure;hist(p_mip_thatHas_atr.*100,[0.5:99.5]);title('% MIP cells expressing ATR');
figure;hist(p_mip_thatHas_tkk.*100,[0.5:99.5]);title('% MIP cells expressing TKK');%%MIP

figure;hist(p_atr_thatHas_fmrf.*100,[0.5:99.5]);title('% ATR cells expressing FMRF');%%ATR
figure;hist(p_atr_thatHas_mip.*100,[0.5:99.5]);title('% ATR cells expressing MIP');
figure;hist(p_atr_thatHas_tkk.*100,[0.5:99.5]);title('% ATR cells expressing TKK');%%ATR

figure;hist(p_tkk_thatHas_fmrf.*100,[0.5:99.5]);title('% TKK cells expressing FMRF');%%TKK
figure;hist(p_tkk_thatHas_mip.*100,[0.5:99.5]);title('% TKK cells expressing MIP');
figure;hist(p_tkk_thatHas_atr.*100,[0.5:99.5]);title('% TKK cells expressing ATR');%%TKK

figure;hist(p_fmrf_thatHas_tkk.*100,[0.5:99.5]);title('% FMRF cells expressing TKK');%%FMRF
figure;hist(p_fmrf_thatHas_mip.*100,[0.5:99.5]);title('% FMRF cells expressing MIP');
figure;hist(p_fmrf_thatHas_atr.*100,[0.5:99.5]);title('% FMRF cells expressing ATR');%%FMRF

close(h);

save FMRF_ATR.mat
%% the averaging script can be used to determine the average % co-expression and standard deviation of every pairwise 
%%co-expression relationship
