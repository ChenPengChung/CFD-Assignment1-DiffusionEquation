# Gnuplot script for grid convergence analysis
set terminal png enhanced size 1920,1080
set output '9.accuracy compare.png'
set title '9.accuracy compare : Difference method- L1 Error vs Grid Spacing'
set xlabel 'log(dx)'
set ylabel 'log(L1 Error)'
set grid
set key left top

# 定義所有擬合函數
f1(x) = 1.92982 * x + -2.54486 
f2(x) = 1.9841 * x + -2.10977
f3(x) = 3.92477 * x + 0.704053
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
    '6.6orderninepoint.dat' using 3:5 pt 6 ps 1.5 lc 'blue' title sprintf('6orderninepoint (slope = %.2f)', 5.24533), \
    f6(x) with lines lw 2 lc 'blue' title sprintf('ninepoint-6order linear (slope = %.2f)', 5.24533)