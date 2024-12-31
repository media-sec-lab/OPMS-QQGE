function [P,mark] = getP_qgm_SI(lambda,element,ary_num,type,e)

mark=[];
if type==1  %for cost
    beta = lambda.*sqrt(pi/2)./element;
elseif type==2   %for standard deviation
    beta = lambda.^(1/4).*element./sqrt(2);
end
if ary_num==3
%     divider_3=normcdf(1.5./beta)-normcdf(-1.5./beta);
%     pChangeP1_3 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta))./divider_3;
%     pChangeM1_3 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta))./divider_3;
% %     pUnChange_3 = (normcdf(0.5./beta)-normcdf(-0.5./beta))./divider_3;
%     pUnChange_3 = 1 - pChangeP1_3 - pChangeM1_3;
%     P=[pChangeM1_3(:),pUnChange_3(:),pChangeP1_3(:)];
%     P=P';
    pChangeP1_3 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta));
    pChangeM1_3 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta));
    pUnChange_3 = (normcdf((0.5-e)./beta)-normcdf((-0.5-e)./beta));
    psum=pChangeP1_3+pChangeM1_3+pUnChange_3;
    pChangeP1_3n = pChangeP1_3./psum;
    pChangeM1_3n = pChangeM1_3./psum;
    pUnChange_3n = pUnChange_3./psum;
    P=[pChangeM1_3n(:),pUnChange_3n(:),pChangeP1_3n(:)];
    P=P';
elseif ary_num==4
    pChangeP2_5 = (normcdf((2.5-e)./beta)-normcdf((1.5-e)./beta));
    pChangeM2_5 = (normcdf((-1.5-e)./beta)-normcdf((-2.5-e)./beta));
    pChangeP1_5 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta));
    pChangeM1_5 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta));
    pUnChange_5 = (normcdf((0.5-e)./beta)-normcdf((-0.5-e)./beta)); 
    mark = pChangeP2_5>=pChangeM2_5;
    pChangePM2 = pChangeP2_5;
    pChangePM2(mark==0) = pChangeM2_5(mark==0);
    psum = pChangeP1_5+pChangeM1_5+pUnChange_5+pChangePM2;    
    pChangeP1_4n = pChangeP1_5./psum;
    pChangeM1_4n = pChangeM1_5./psum;
    pUnChange_4n = pUnChange_5./psum;
    pChangePM2_4n = pChangePM2./psum;
    P=[pChangeM1_4n(:),pUnChange_4n(:),pChangeP1_4n(:),pChangePM2_4n(:)];
    P=P';
elseif ary_num==5
%     divider_5=normcdf(2.5./beta)-normcdf(-2.5./beta);
%     pChangeP2_5 = (normcdf((2.5-e)./beta)-normcdf((1.5-e)./beta))./divider_5;
%     pChangeM2_5 = (normcdf((-1.5-e)./beta)-normcdf((-2.5-e)./beta))./divider_5;
%     pChangeP1_5 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta))./divider_5;
%     pChangeM1_5 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta))./divider_5;
% %     pUnChange_5=(normcdf(0.5./beta)-normcdf(-0.5./beta))./(normcdf(2.5./beta)-normcdf(-2.5./beta));    
%     pUnChange_5 = 1 - pChangeP2_5 - pChangeM2_5 - pChangeP1_5 - pChangeM1_5;  
%     P=[pChangeM2_5(:),pChangeM1_5(:),pUnChange_5(:),pChangeP1_5(:),pChangeP2_5(:)];
%     P=P';
    pChangeP2_5 = (normcdf((2.5-e)./beta)-normcdf((1.5-e)./beta));
    pChangeM2_5 = (normcdf((-1.5-e)./beta)-normcdf((-2.5-e)./beta));
    pChangeP1_5 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta));
    pChangeM1_5 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta));
    pUnChange_5 = (normcdf((0.5-e)./beta)-normcdf((-0.5-e)./beta)); 
    psum = pChangeP1_5+pChangeM1_5+pUnChange_5+pChangeP2_5+pChangeM2_5;    
    pChangeP1_5n = pChangeP1_5./psum;
    pChangeM1_5n = pChangeM1_5./psum;
    pUnChange_5n = pUnChange_5./psum;
    pChangeP2_5n = pChangeP2_5./psum;
    pChangeM2_5n = pChangeM2_5./psum;
    P=[pChangeM2_5n(:),pChangeM1_5n(:),pUnChange_5n(:),pChangeP1_5n(:),pChangeP2_5n(:)];
    P=P';
elseif ary_num==7
%     divider_7=normcdf(3.5./beta)-normcdf(-3.5./beta);
%     pChangeP3_7 = (normcdf((3.5-e)./beta)-normcdf((2.5-e)./beta))./divider_7;
%     pChangeM3_7 = (normcdf((-2.5-e)./beta)-normcdf((-3.5-e)./beta))./divider_7;
%     pChangeP2_7 = (normcdf((2.5-e)./beta)-normcdf((1.5-e)./beta))./divider_7;
%     pChangeM2_7 = (normcdf((-1.5-e)./beta)-normcdf((-2.5-e)./beta))./divider_7;
%     pChangeP1_7 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta))./divider_7;
%     pChangeM1_7 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta))./divider_7;
% %     pUnChange_5=(normcdf(0.5./beta)-normcdf(-0.5./beta))./(normcdf(2.5./beta)-normcdf(-2.5./beta));    
%     pUnChange_7 = 1 - pChangeP3_7 - pChangeM3_7 - pChangeP2_7 - pChangeM2_7 - pChangeP1_7 - pChangeM1_7;  
%     P=[pChangeM3_7(:),pChangeM2_7(:),pChangeM1_7(:),pUnChange_7(:),pChangeP1_7(:),pChangeP2_7(:),pChangeP3_7(:)];
%     P=P';
    pChangeP3_7 = (normcdf((3.5-e)./beta)-normcdf((2.5-e)./beta));
    pChangeM3_7 = (normcdf((-2.5-e)./beta)-normcdf((-3.5-e)./beta));
    pChangeP2_7 = (normcdf((2.5-e)./beta)-normcdf((1.5-e)./beta));
    pChangeM2_7 = (normcdf((-1.5-e)./beta)-normcdf((-2.5-e)./beta));
    pChangeP1_7 = (normcdf((1.5-e)./beta)-normcdf((0.5-e)./beta));
    pChangeM1_7 = (normcdf((-0.5-e)./beta)-normcdf((-1.5-e)./beta));
    pUnChange_7 = (normcdf((0.5-e)./beta)-normcdf((-0.5-e)./beta));    
    psum = pChangeP3_7+pChangeM3_7+pChangeP2_7+pChangeM2_7+pChangeP1_7+pChangeM1_7+pUnChange_7;    
    pChangeP1_7n = pChangeP1_7./psum;
    pChangeM1_7n = pChangeM1_7./psum;
    pUnChange_7n = pUnChange_7./psum;
    pChangeP2_7n = pChangeP2_7./psum;
    pChangeM2_7n = pChangeM2_7./psum;  
    pChangeP3_7n = pChangeP3_7./psum;
    pChangeM3_7n = pChangeM3_7./psum;     
    P=[pChangeM3_7n(:),pChangeM2_7n(:),pChangeM1_7n(:),pUnChange_7n(:),pChangeP1_7n(:),pChangeP2_7n(:),pChangeP3_7n(:)];
    P=P';
end

end

