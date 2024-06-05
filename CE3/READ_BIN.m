
clear all ; 
close all ; 
clc ; 

[ y_rst , r_rst , u_rst , d_rst , anyRef_rst, sample_rst ] = read_binary_file('logQ4cm.bin') ; 
[ y_dat1 , r_dat1 , u_dat1 , d_dat1 , anyRef_dat1 , sample_dat1 ] = read_binary_file('log4cmdatadriven.bin') ; 
[ y_dat2, r_dat2 , u_dat2 , d_dat2 , anyRef_dat2 , sample_dat2 ] = read_binary_file('log8cmdatadriven.bin') ; 
[ y_dat3, r_dat3 , u_dat3 , d_dat3 , anyRef_dat3 , sample_dat3 ] = read_binary_file('log14cmdatadriven.bin') ; 


%%

figure(1) 
plot( r_rst(:,1) )
hold on
plot( y_rst(: ,1) )
hold off 
title( "output tracking a reference using an RST controller" )
legend('reference' , 'output' )

figure(2) 
plot( u_rst(: ,1) )
title( "input computed using our RST controller" )

%%

figure(3) 
plot( r_dat1(:,1) )
hold on
plot( y_dat1(: ,1) )
hold off 
title("output tracking a reference using an datadriven controller 1")
legend('reference','output')

figure(4) 
plot( u_dat1(: ,1) )
title("input computed using our datadriven controller 1")

%%

figure(5) 
plot( r_dat2(:,1) )
hold on
plot( y_dat2(: ,1) )
hold off 
title("output tracking a reference using an datadriven controller 2")
legend('reference','output')

figure(6) 
plot( u_dat2(: ,1) )
title("input computed using our datadriven controller 2")

%%

figure(7) 
plot( r_dat3(:,1) )
hold on
plot( y_dat3(: ,1) )
hold off 
title("output tracking a reference using an datadriven controller 3")
legend('reference','output')

figure(8) 
plot( u_dat3(: ,1) )
title("input computed using our datadriven controller 3")
