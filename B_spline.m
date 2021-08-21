function [curve,curve1,curve2,B_3]=B_spline(sample,n,d,u,resolution,type)
    if(d<1)
        fprintf("输入次数小于1");
        return
    end
    
    switch type
        case 'uniform'
             curve_size = u(n+d)/resolution;
             interval = (u(2)-u(1))/resolution;
        case 'quasi-uniform'
             curve_size = (u(d+1)-u(d))/resolution*(n-d+1);
             interval = (u(d+1)-u(d))/resolution;
    end

    d_now = 2;
    %求各段函数的参数a、b
    B_2(1) = struct;
    B_2(1).section1 = struct;
    B_2(1).section2 = struct;
    B_2(1).section1.a=0;
    B_2(1).section1.b=0;
    B_2(1).section2.a=0;
    B_2(1).section2.b=0;
    for k=1:(n+d-2)
        %各个分段函数的次数，section1(0)为常数项，section1(1)为一次项系数
        %section1的定义域为(u(k),u(k+1))，section2的定义域为(u(k+1),u(k+2))
        B_2(k).section1.a=0;
        B_2(k).section1.b=0;
        B_2(k).section2.a=0;
        B_2(k).section2.b=0;

        if (u(k+d_now-1)-u(k))~=0
            B_2(k).section1.b = -u(k)/(u(k+d_now-1)-u(k));
            B_2(k).section1.a = 1/(u(k+d_now-1)-u(k));
        end

        if (u(k+d_now)-u(k+1))~=0
            B_2(k).section2.b = u(k+d_now)/(u(k+d_now)-u(k+1));
            B_2(k).section2.a = -1/(u(k+d_now)-u(k+1));
        end
    end
    
    switch type
        case 'uniform'
            %拟合1次曲线
            curve1 = zeros(n+d-2,curve_size);
            for k=1:(n+d-2)
                for i=1:(2*interval)
                    if i<=interval
                        curve1(k,i+(k-1)*interval) = curve1(k,i+(k-1)*interval) + B_2(k).section1.a * (u(k)+i*resolution) + B_2(k).section1.b;
                    else
                        curve1(k,i+(k-1)*interval) = curve1(k,i+(k-1)*interval) + B_2(k).section2.a * (u(k)+i*resolution) + B_2(k).section2.b;
                    end
                end
            end
            if(d==2)
                curve = zeros(1,curve_size);
                for k = 1:n
                    curve = curve + sample(k) * curve1(k,:);
                end
                return
            end
        case 'quasi-uniform'
             curve1 = zeros(n-d+2,curve_size);
             for i = 1:interval
                 curve1(1,i) = B_2(d-1).section2.a * (u(d)+i*resolution) + B_2(d-1).section2.b;
             end
             for k = 1:(n-d)
                 for i=1:(2*interval)
                    if i<=interval
                        curve1(k+1,i+(k-1)*interval) =  B_2(d-1+k).section1.a * (u(d-1+k)+i*resolution) + B_2(d-1+k).section1.b;
                    else
                        curve1(k+1,i+(k-1)*interval) =  B_2(d-1+k).section2.a * (u(d-1+k)+i*resolution) + B_2(d-1+k).section2.b;
                    end
                end
             end
             for i = 1:interval
                 curve1(n+d-2,i+k*interval) = B_2(d+k).section1.a * (u(d-1+k)+i*resolution) + B_2(d-1+k).section1.b;
             end
             
             if(d==2)
                curve = zeros(1,curve_size);
                curve2 = zeros(1,1);
                B_3 = 0;
                for k = 1:n
                    curve = curve + sample(k) * curve1(k,:);
                end
                return
            end
    end

    % a*u^2+b*u+c
    d_now = 3;
    B_3(1) = struct;
    B_3(1).section1 = struct;
    B_3(1).section2 = struct;
    B_3(1).section3 = struct;
    B_3(1).section1.a=0;
    B_3(1).section1.b=0;
    B_3(1).section1.c=0;
    B_3(1).section2.a=0;
    B_3(1).section2.b=0;
    B_3(1).section2.c=0;
    B_3(1).section3.a=0;
    B_3(1).section3.b=0;
    B_3(1).section3.c=0;

    for k=1:(n+d-3)
        B_3(k).section1.a=0;
        B_3(k).section1.b=0;
        B_3(k).section1.c=0;
        B_3(k).section2.a=0;
        B_3(k).section2.b=0;
        B_3(k).section2.c=0;
        B_3(k).section3.a=0;
        B_3(k).section3.b=0;
        B_3(k).section3.c=0;

        if (u(k+d_now-1)-u(k))~=0
            B_3(k).section1.a = B_2(k).section1.a/(u(k+d-1)-u(k));
            B_3(k).section1.b = (B_2(k).section1.b-u(k)*B_2(k).section1.a)/(u(k+d-1)-u(k));
            B_3(k).section1.c = -u(k)*B_2(k).section1.b/(u(k+d-1)-u(k));

            B_3(k).section2.a = B_3(k).section2.a + B_2(k).section2.a/(u(k+d-1)-u(k));
            B_3(k).section2.b = B_3(k).section2.b + (B_2(k).section2.b-B_2(k).section2.a*u(k))  /(u(k+d-1)-u(k));
            B_3(k).section2.c = B_3(k).section2.c - u(k)*B_2(k).section2.b/(u(k+d-1)-u(k));
        end

        if (u(k+d_now)-u(k+1))~=0
            B_3(k).section2.a = B_3(k).section2.a + (-B_2(k+1).section1.a)/(u(k+d_now)-u(k+1));
            B_3(k).section2.b = B_3(k).section2.b + (u(k+d_now)*B_2(k+1).section1.a-B_2(k+1).section1.b)/(u(k+d_now)-u(k+1));
            B_3(k).section2.c = B_3(k).section2.c + (u(k+d_now)*B_2(k+1).section1.b)/(u(k+d_now)-u(k+1));

            B_3(k).section3.a = (-B_2(k+1).section2.a)/(u(k+d_now)-u(k+1));
            B_3(k).section3.b = (u(k+d_now)*B_2(k+1).section2.a-B_2(k+1).section2.b)/(u(k+d_now)-u(k+1));
            B_3(k).section3.c = u(k+d_now)*B_2(k+1).section2.b/(u(k+d_now)-u(k+1));
        end
    end

    switch type
        case 'uniform'
            curve2 = zeros(n+d-3,curve_size);
            %拟合2次曲线
            for k=1:(n+d-3)
                for i=1:(3*interval)
                    if i<=interval
                        curve2(k,i+(k-1)*interval) = curve2(k,i+(k-1)*interval) + B_3(k).section1.a * (u(k)+i*resolution)^2 + B_3(k).section1.b*(u(k)+i*resolution) + B_3(k).section1.c;
                    elseif i<=2*interval && i>interval
                        curve2(k,i+(k-1)*interval) = curve2(k,i+(k-1)*interval) + B_3(k).section2.a * (u(k)+i*resolution)^2 + B_3(k).section2.b*(u(k)+i*resolution) + B_3(k).section2.c;
                    else
                        curve2(k,i+(k-1)*interval) = curve2(k,i+(k-1)*interval) + B_3(k).section3.a * (u(k)+i*resolution)^2 + B_3(k).section3.b*(u(k)+i*resolution) + B_3(k).section3.c;
                    end
                end
            end

            if(d==3)
                curve = zeros(1,curve_size);
                for k = 1:n
                    curve = curve + sample(k) * curve2(k,:);
                end
                return
            end
            
        case 'quasi-uniform'
            curve2 = zeros(n-d+3,curve_size);
            %拟合2次曲线
            for i=1:interval
                curve2(1,i) = B_3(d-2).section3.a * (u(d)+i*resolution)^2 + B_3(d-2).section3.b*(u(d)+i*resolution) + B_3(d-2).section3.c;
                
                curve2(2,i) = B_3(d-1).section2.a * (u(d)+i*resolution)^2 + B_3(d-1).section2.b*(u(d)+i*resolution) + B_3(d-1).section2.c;
                curve2(2,interval+i) = B_3(d-1).section3.a * (u(d+1)+i*resolution)^2 + B_3(d-1).section3.b*(u(d+1)+i*resolution) + B_3(d-1).section3.c;
                
                curve2(n-d+2,(n-d-1)*interval+i) = B_3(n-1).section1.a * (u(n-1)+i*resolution)^2 + B_3(n-1).section1.b*(u(n-1)+i*resolution) + B_3(n-1).section1.c;
                curve2(n-d+2,(n-d)*interval+i) = B_3(n-1).section2.a * (u(n)+i*resolution)^2 + B_3(n-1).section2.b*(u(n)+i*resolution) + B_3(n-1).section2.c;
                
                curve2(n-d+3,(n-d)*interval+i) = B_3(n).section1.a * (u(n)+i*resolution)^2 + B_3(n).section1.b*(u(n)+i*resolution) + B_3(n).section1.c;
            end
            
            for k=1:(n-d-1)
                for i=1:(3*interval)
                    if i<=interval
                        curve2(k+2,i+(k-1)*interval) = B_3(d-1+k).section1.a * (u(d-1+k)+i*resolution)^2 + B_3(d-1+k).section1.b*(u(d-1+k)+i*resolution) + B_3(d-1+k).section1.c;
                    elseif i<=2*interval && i>interval
                        curve2(k+2,i+(k-1)*interval) = B_3(d-1+k).section2.a * (u(d-1+k)+i*resolution)^2 + B_3(d-1+k).section2.b*(u(d-1+k)+i*resolution) + B_3(d-1+k).section2.c;
                    else
                        curve2(k+2,i+(k-1)*interval) = B_3(d-1+k).section3.a * (u(d-1+k)+i*resolution)^2 + B_3(d-1+k).section3.b*(u(d-1+k)+i*resolution) + B_3(d-1+k).section3.c;
                    end
                end
            end
            
            if(d==3)
                curve = zeros(1,curve_size);
                for k = 1:n
                    curve = curve + sample(k) * curve2(k,:);
                end
                return
            end
            
    end

            

    

end

