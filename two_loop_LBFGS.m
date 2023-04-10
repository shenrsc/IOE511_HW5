function [d] =  two_loop_LBFGS(m,s_all,y_all,f,g,h0k)
    q = g;
    size_s = size(s_all);
    if(size_s(2)>0)
        rou = zeros(1,size_s(2));
        alpha = zeros(1,size_s(2));
    end

    for i=size_s(2):-1:1
        rou(i) = 1/(s_all(:,i).'*y_all(:,i));
        alpha(i) = rou(i)*s_all(:,i).'*q;
        q = q-alpha(i)*y_all(:,i);
    end

    r = h0k*q;
    for i=1:1:size_s(2)
        beta = rou(i)*y_all(:,i).'*r;
        r = r+s_all(:,i)*(alpha(i)-beta);
    end
    d = -r;
end