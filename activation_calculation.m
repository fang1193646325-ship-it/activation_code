%-----------------------主函数----------------------------
function main()
    
    %activation_cal_test();
    activation_cal();
    

end

%-----------------------子函数-----------------------------


%（1）计算活度【未完成】
function activation_cal()
    
    time_step = 300;%设置时间步长为300s(有待商榷)
    
    %中子通量计算【输入变量】
    phi = 1*1e12;
    
    %所有核素反应截面数据(单位：cm^2)【输入变量】
    sigma_Ni58 = 1e-25;
    sigma_Ni59 = 1e-25;
    sigma_Co59 = 1e-26;
    

    %所有核素衰变常数
    lambda_Ni59 = log(2)/(76002*365*24*60*60);
    lambda_Co60 = log(2)/(5.2712*365*24*60*60);

    %所有核素数量大小【输入变量】
    number_Ni58 = 1e20;
    number_Ni59 = 0;
    number_Ni60 = 1e20;
    number_Co59 = 1e20;
    number_Co60 = 0;

    %初始化A矩阵【待优化】
    A = [-sigma_Ni58*phi,0,0,0,0;
         sigma_Ni58*phi,-(sigma_Ni59*phi+lambda_Ni59),0,0,0;
         0,sigma_Ni59*phi,0,0,lambda_Co60;
         0,lambda_Ni59,0,-sigma_Co59*phi,0;
         0,0,0,sigma_Co59*phi,-lambda_Co60
         ];

    %初始核素条件【待优化】
    N0 = [number_Ni58;
          number_Ni59;
          number_Ni60;
          number_Co59;
          number_Co60];
    
    Nt = [];%列表示每个步长，行表示每种核素
    T = [];
    for t=0:time_step:60*60*24*365*1%计算15a内的测试数据
        
        %N = taylor_calculation(A,N0,t);
        N = Chebyshev_calculation(A,N0,t)
        T = [T,t];
        Nt = [Nt,N];
     
    end
    
    Nt = [T;Nt];
    
    Drawing_Picture(Nt,T);
    
    Nt = Nt.';

    fid = fopen('test_data.txt','w');
    fprintf('时间 Ni58 Ni59 Ni60 Co59 Co60');
    for i = 1:size(Nt(:,1))
        fprintf(fid,'%d %d %d %d %d %d\n',Nt(i,:));
    end
    fclose(fid);

end
%（辅助功能）计算活化【测试函数】
function activation_cal_test()
    time_step = 0.1;%设置时间步长为0.1s(有待商榷)
    
    %初始化A矩阵(待完善)
    A = [1,2;
         3,4];

    %初始核素条件
    N0 = [1e20;
          1e20];
    
    test = [];
    T = [];
    clc;
    for t=0:time_step:1%计算1s内的测试数据
        B = Chebyshev_calculation(A,N0,t);
        T = [T;t;t];
        test = [test;B];
    end
    
    test = [T,test];
    xlswrite('test_data.xlsx',test);
end

%（辅助功能）画图函数
function Drawing_Picture(Nt,T)
     %画图
        %设置变量
        x = T;
        y_Ni58 = Nt(2,:);
        y_Ni59 = Nt(3,:);
        y_Ni60 = Nt(4,:);
        y_Co59 = Nt(5,:);
        y_Co60 = Nt(6,:);
        
        figure;

        subplot(2,1,1)
        %绘制折线图1，设置线型、颜色和标记
        plot(x,y_Ni58,'-r',x,y_Ni60,'-g',x,y_Co59,'-b');
        %右上角标注
        legend('Ni58','Ni60','Co59');
        %x轴坐标描述
        xlabel('时间(s)');
        %y轴坐标描述
        ylabel('核素数量');
        
        subplot(2,1,2);
        %绘制折线图2，设置线型、颜色和标记
        plot(x,y_Ni59,'-c',x,y_Co60,'-y');
        %右上角标注
        legend('Ni59','Co60');
        %x轴坐标描述
        xlabel('时间(s)');
        %y轴坐标描述
        ylabel('核素数量');

end

%（1.1）计算当前的中子通量,更新A矩阵【未完成】
function A = flux_calculation(A)
    
    %待完善

end 

%（1.2）泰勒展开法计算下一时刻的核素数量矩阵
function Nt = taylor_calculation(A,N0,t)%A为方程组的系数矩阵、N0为初始核素矩阵，t为某一时刻
    
     %确定泰勒展开到第几项级数，n代表泰勒展开到第n项
    %n = fix(min(norm(A,inf),norm(A,1))*t);
    n = 10;
    
    %为了减少重复计算，设置三个临时变量存储A^n、t^n、k!
    %考虑精度问题，需要用double变量去存储所有数据
    temp_A = A;%A矩阵的幂次方
    temp_t = t;%时间t的幂次方
    temp_k = 1.0;%k！阶乘

    %初始化e^At矩阵
    Nt = eye(size(A));
        
    %从泰勒展开第一项到第n项相加
    for k=1:n
        temp_k = temp_k*k;
        Nt = Nt + temp_A*temp_t/temp_k;%每一项累加
        temp_A = temp_A*A;
        temp_t = temp_t*t;
    end    
    
    Nt = Nt * N0;

    return;
end

%（1.3）切比雪夫法计算下一时刻的核素数量矩阵
function Nt = Chebyshev_calculation(A,N0,t)%A为方程组的系数矩阵、N0为初始核素矩阵，t为某一时刻
    
    %a0复数
    alpha0 = 1.8321743782540412751e-14;

    %alpha的实部参数
    alpha_Real_part = [
         -7.1542880635890672853e-5,9.4390253107361688779e-3,-3.7636003878226968717e-1,-2.3498232091082701191e1,4.6933274488831293047e1,-2.7875161940145646468e1,4.8071120988325088907e0
    ];

    %alpha的虚部参数
    alpha_Imaginary_part = [
       1.4361043349541300111e-4,-1.7184791958483017511e-2,3.3518347029450104214e-1,-5.8083591297142074004e0,4.5643649768827760791e1,-1.0214733999056451434e2,-1.3209793837428723881e0
    ];

    %theta的实部参数
    theta_Real_part = [
        -8.8977731864688888199e0,-3.7032750494234480603e0,-0.2087586382501301251e0,3.9933697105785685194e0,5.0893450605806245066e0,5.6231425727459771248e0,2.2697838292311127097e0
    ];

    %theta的虚部参数
    theta_Imaginary_part = [
        1.6630982619902085304e1,1.3656371871483268171e1,1.0991260561901260913e1,6.0048316422350373178e0,3.5888240290270065102e0,1.1940690463439669766e0,8.4617379730402214019e0
    ];
    
    %构建alpha复数
    alpha = complex(alpha_Real_part,alpha_Imaginary_part);
    %构建theta复数
    theta = complex(theta_Real_part,theta_Imaginary_part);
   
    I = eye(size(A));%A矩阵规模的单位矩阵
    sum = 0*N0;
    A = A*t
    
    
    for i=1:7
        sum = sum + (A - theta(i)*I) \ (alpha(i)*N0);%累加
    end
    
    Nt = alpha0*N0 + 2*real(sum);%a0-2*Re（求和的实部）

    return 
end
 
