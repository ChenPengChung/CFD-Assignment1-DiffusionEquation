# Gnuplot script for grid convergence analysis
set terminal png enhanced size 1920,1080
set output '8.9 points comapre.png'
set title '8.9 points comapre : Difference method- L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top

# 定義所有擬合函數
f4(x) = 3.90919 * x + 0.532055
f5(x) = 3.97604 * x + 0.234779
f6(x) = 5.24533 * x + -6.10883

# 使用單一個 plot 指令繪製所有數據和擬合線
plot \
    '4.ninepoint-2order.dat' using 3:5 pt 8 ps 1.5 lc 'red' title sprintf('ninepoint-2order (slope = %.2f)', 3.90919), \
    f4(x) with lines lw 2 lc 'red' title sprintf('ninepoint-2order linear (slope = %.2f)', 3.90919), \
    \
    '5.ninepoint-4order.dat' using 3:5 pt 5 ps 1.5 lc 'purple' title sprintf('ninepoint-4order (slope = %.2f)', 3.97604), \
    f5(x) with lines lw 2 lc 'purple' title sprintf('ninepoint-4order linear (slope = %.2f)', 3.97604), \
    \
    '6.6orderninepoint.dat' using 3:5 pt 6 ps 1.5 lc 'blue' title sprintf('6orderninepoint (slope = %.2f)', 5.24533), \
    f6(x) with lines lw 2 lc 'blue' title sprintf('ninepoint-6order linear (slope = %.2f)', 5.24533)