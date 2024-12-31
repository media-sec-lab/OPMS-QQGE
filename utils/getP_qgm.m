function [P] = getP_qgm(lambda,element,ary_num,type)

if type==1  %for cost
    beta = lambda.*sqrt(pi/2)./element;
elseif type==2   %for standard deviation
    beta = lambda.^(1/4).*element./sqrt(2);
end
if ary_num==3
    divider_3=normcdf(1.5./beta)-normcdf(-1.5./beta);
    pChangeP1_3 = (normcdf(1.5./beta)-normcdf(0.5./beta))./divider_3;
    pChangeM1_3 = (normcdf(-0.5./beta)-normcdf(-1.5./beta))./divider_3;
%     pUnChange_3 = (normcdf(0.5./beta)-normcdf(-0.5./beta))./divider_3;
    pUnChange_3 = 1 - pChangeP1_3 - pChangeM1_3;
    P=[pChangeM1_3(:),pUnChange_3(:),pChangeP1_3(:)];
    P=P';
elseif ary_num==5
    divider_5=normcdf(2.5./beta)-normcdf(-2.5./beta);
    pChangeP2_5 = (normcdf(2.5./beta)-normcdf(1.5./beta))./divider_5;
    pChangeM2_5 = (normcdf(-1.5./beta)-normcdf(-2.5./beta))./divider_5;
    pChangeP1_5 = (normcdf(1.5./beta)-normcdf(0.5./beta))./divider_5;
    pChangeM1_5 = (normcdf(-0.5./beta)-normcdf(-1.5./beta))./divider_5;
%     pUnChange_5=(normcdf(0.5./beta)-normcdf(-0.5./beta))./(normcdf(2.5./beta)-normcdf(-2.5./beta));    
    pUnChange_5 = 1 - pChangeP2_5 - pChangeM2_5 - pChangeP1_5 - pChangeM1_5;  
    P=[pChangeM2_5(:),pChangeM1_5(:),pUnChange_5(:),pChangeP1_5(:),pChangeP2_5(:)];
    P=P';
elseif ary_num==7
    divider_7=normcdf(3.5./beta)-normcdf(-3.5./beta);
    pChangeP3_7 = (normcdf(3.5./beta)-normcdf(2.5./beta))./divider_7;
    pChangeM3_7 = (normcdf(-2.5./beta)-normcdf(-3.5./beta))./divider_7;
    pChangeP2_7 = (normcdf(2.5./beta)-normcdf(1.5./beta))./divider_7;
    pChangeM2_7 = (normcdf(-1.5./beta)-normcdf(-2.5./beta))./divider_7;
    pChangeP1_7 = (normcdf(1.5./beta)-normcdf(0.5./beta))./divider_7;
    pChangeM1_7 = (normcdf(-0.5./beta)-normcdf(-1.5./beta))./divider_7;
%     pUnChange_5=(normcdf(0.5./beta)-normcdf(-0.5./beta))./(normcdf(2.5./beta)-normcdf(-2.5./beta));    
    pUnChange_7 = 1 - pChangeP3_7 - pChangeM3_7 - pChangeP2_7 - pChangeM2_7 - pChangeP1_7 - pChangeM1_7;  
    P=[pChangeM3_7(:),pChangeM2_7(:),pChangeM1_7(:),pUnChange_7(:),pChangeP1_7(:),pChangeP2_7(:),pChangeP3_7(:)];
    P=P';
end

end

