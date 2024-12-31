function [lambda, p, mark] = calc_lambda_curve_forward_lt1_qgm_SI(element, message_length, n, ary_num,type, e)
%% Monotonic increasing    

% m3 = double(message_length + 1);
m3 = double(message_length - 1);
l1=0;
l3=0;
%         m3=double(n);
%         m3=4.154884097986136e+05;
iterations = 0;
syms x;
sta=0.1;    %¿É¸ù¾ÝÊµ¼Ê½øÐÐÉè¶¨


% while m3 > message_length
while m3 < message_length    %ÕÒµ½ÓÒ¶Ëµã
    l1=l3;
    m1=m3;
    l3 = l3 * 2;            
    if l3==0
        l3=sta;
    end
    [p,mark] = getP_qgm_SI(l3,element,ary_num,type,e);
    m3 = getHt_qgm(p);
    iterations = iterations + 1;
    if (iterations > 300)
        lambda = l3;
        return;
    end
end      
                             
lambda=sta; 
lastm=0;
        
if l3==sta
    while (double(abs(lastm/n-message_length/n)) > 1/1000.0 ) && (iterations<300)
        lambda = l1+(l3-l1)/2; 
        [p,mark] = getP_qgm_SI(lambda,element,ary_num,type,e);
        m2 = getHt_qgm(p);                       
%         if m2 < message_length
        if m2 > message_length
            l3 = lambda;
            m3 = m2;
            lastm = m3;
        else
            l1 = lambda;
            m1 = m2;
            lastm=m1;
        end
        iterations = iterations + 1;
    end         
else
    k=(m3-m1)/(l3-l1);   %¶þ·Ö·¨ÖÐµã
    l2=(message_length-m1)/k+l1;
    [p,mark] = getP_qgm_SI(l2,element,ary_num,type,e);
    m2 = getHt_qgm(p);
    iterations = iterations + 1;
    lambda=l2;
    while  (double(abs(lastm/n-message_length/n)) > 1/1000.0 ) && (iterations<300)
        z=polyfit([l1 l2 l3],[m1,m2,m3],2);   %¶þ´ÎÇúÏßÄâºÏ
        xz=double(solve(z(1)*x*x+z(2)*x+z(3)-message_length));
        if xz(1)>=l1 && xz(1)<=l3
            lambda=xz(1);
        else
            lambda=xz(2);
        end
        [p,mark] = getP_qgm_SI(lambda,element,ary_num,type,e);
        lastm = getHt_qgm(p);
        iterations = iterations + 1;  
        if double(abs(lastm/n-message_length/n)) < 1/1000.0
            break;
        end     
        if m2>lastm && m2>message_length
%             l1=l2;  m1=m2;  l2=lambda;  m2=lastm;
            l3=l2;  m3=m2;  l2=lambda;  m2=lastm;            
        end
        if m2<lastm && m2<message_length
%             l3=l2;  m3=m2;  l2=lambda;  m2=lastm;
            l1=l2;  m1=m2;  l2=lambda;  m2=lastm;
        end
        if lastm>m2 && m2>message_length
%             l1=lambda;  m1=lastm;
            l3=lambda;  m3=lastm;            
        end
        if message_length>m2 && m2>lastm
%             l3=lambda;  m3=lastm;
            l1=lambda;  m1=lastm;
        end
    end
end   

end
