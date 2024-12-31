function [lambda, p, mark] = calc_lambda_qgm_SI(element, message_length, n, ary_num, type, e)
%% Monotonic increasing
        l3 = 1e+3;
%         m3 = double(message_length + 1);
        m3 = double(message_length - 1);
        iterations = 0;
        p = zeros(size(element));
%         while m3 > message_length
        while m3 < message_length    %�ҵ��Ҷ˵�
            l3 = l3 * 2;
%             p = getP_qgm(l3,element,ary_num,type);
            [p,mark] = getP_qgm_SI(l3,element,ary_num,type,e);
            m3 = getHt_qgm(p);
            iterations = iterations + 1;
            if (iterations > 300)
                lambda = l3;
                return;
            end
        end        
        
        l1 = 0; 
        m1 = double(n);        
        lambda = 0;
        
        alpha = double(message_length)/n;
%         disp(message_length/n);
        % limit search to 30 iterations
        % and require that relative payload embedded is roughly within 1/1000 of the required relative payload   
        lastm=0;
        while  (double(abs(lastm/n-message_length/n)) > 1/1000.0 ) && (iterations<300)
%         while  ((m1 - m3) / n > alpha / n * 1e-2)
            lambda = l1+(l3-l1)/2; 
%             p = getP_qgm(lambda,element,ary_num,type);
            [p,mark] = getP_qgm_SI(lambda,element,ary_num,type,e);
            m2 = getHt_qgm(p);
%     		if m2 < message_length
            if m2 > message_length
    			l3 = lambda;
    			m3 = m2;
                lastm = m3;
            else
    			l1 = lambda;
    			m1 = m2;
                lastm=m1;
            end
%             disp(lastm/n)
    		iterations = iterations + 1;
%             disp(iterations);
        end
end