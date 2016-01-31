set terminal pdf
set output 'all_fe.pdf'

#fe
fi='tmp2.csv'
set xlabel 'Temperature (°C)'
set datafile separator ','
set title 'Correlation Length by Temperature'
set ylabel '\xi / r_{mean} (nm)'
plot fi u 2:3 w lp t 'x direction', fi u 2:9 w lp t 'y direction', fi u 2:15 w lp t 'mean grain radius'

set title 'Interface Layer Thickness by Temperature'
set ylabel 'ILT (nm)'
plot fi u 2:4 w lp t 'x direction', fi u 2:10 w lp t 'y direction'

set title 'Roughness exponent by Temperature'
set ylabel 'RE'
plot fi u 2:($5/2) w lp t 'x direction', fi u 2:($11/2) w lp t 'y direction'

#fe_pt
set output 'all_fept.pdf'
fi='tmppt.csv'
set title 'Correlation Length by Temperature'
set ylabel '\xi / r_{mean} (nm)'
plot fi u 2:3 w lp t 'x direction', fi u 2:9 w lp t 'y direction', fi u 2:15 w lp t 'mean grain radius'

set title 'Interface Layer Thickness by Temperature'
set ylabel 'ILT (nm)'
plot fi u 2:4 w lp t 'x direction', fi u 2:10 w lp t 'y direction'

set title 'Roughness exponent by Temperature'
set ylabel 'RE'
plot fi u 2:($5/2) w lp t 'x direction', fi u 2:($11/2) w lp t 'y direction'