# Gnuplot script for grid convergence analysis
set terminal png enhanced size 1920,1080
set output '7.accuracy for 2-order,4-order,9points.png'
set title '7.accuracy for 2-order,4-order,9points- L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top

# 定義所有擬合函數
f1(x) = 1.92982 * x + -2.54486
f2(x) = 1.9841 * x + -2.10977
f3(x) = 3.92477 * x + 0.704053
f4(x) = 3.90919 * x + 0.532055
f5(x) = 3.97604 * x + 0.234779
f6(x) = 5.24533 * x + -6.10883

# 使用單一個 plot 指令繪製所有數據和擬合線
plot \
    '1.2orderCD_FDM.dat' using 3:5 pt 7 ps 1.5 lc 'green' title sprintf('2orderCD_FDM (slope = %.2f)', 1.92982), \
    f1(x) with lines lw 2 lc 'green' title sprintf('2orderCD_FDM linear (slope = %.2f)', 1.92982), \
    \
    '2.2orderCD_FVM.dat' using 3:5 pt 14 ps 1.5 lc 'orange' title sprintf('2orderCD_FVM (slope = %.2f)', 1.9841), \
    f2(x) with lines lw 2 lc 'orange' title sprintf('2orderCD_FVM linear (slope = %.2f)', 1.9841), \
    \
    '3.4orderCD-2order.dat' using 3:5 pt 13 ps 1.5 lc 'brown' title sprintf('4orderCD-2order (slope = %.2f)', 3.92477), \
    f3(x) with lines lw 2 lc 'brown' title sprintf('4orderCD-2order linear (slope = %.2f)', 3.92477), \
    \
    '4.ninepoint-2order.dat' using 3:5 pt 8 ps 1.5 lc 'red' title sprintf('ninepoint-2order (slope = %.2f)', 3.90919), \
    f4(x) with lines lw 2 lc 'red' title sprintf('ninepoint-2order linear (slope = %.2f)', 3.90919), \
    \
    '5.ninepoint-4order.dat' using 3:5 pt 5 ps 1.5 lc 'purple' title sprintf('ninepoint-4order (slope = %.2f)', 3.97604), \
    f5(x) with lines lw 2 lc 'purple' title sprintf('ninepoint-4order linear (slope = %.2f)', 3.97604), \
    \
    '6.6orderninepoint.dat' using 3:5 pt 6 ps 1.5 lc 'blue' title sprintf('6orderninepoint (slope = %.2f)', 5.24533), \
    f6(x) with lines lw 2 lc 'blue' title sprintf('ninepoint-6order linear (slope = %.2f)', 5.24533)