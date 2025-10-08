set terminal png enhanced size 800,600
set output '5.ninepoint-4order.png'
set title '5.ninepoint-4order : L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top
f(x) = 3.97604 * x + 0.234779
plot 'ninepoint-4order.dat' using 3:5  pt 5 ps 1.5 lc 'purple' title sprintf('ninepoint-4order (slope = %.2f)', 3.97604), \
     f(x) with lines lw 2 lc 'purple' title sprintf('ninepoint-4order linear (slope = %.2f)', 3.97604), \
    
