set terminal png size 700,800 font "Verdana,32"
set output outputFileName

set xlabel "t [s]"
set ylabel "E_{K-tot} [J]"
#set title "Hydrostatický tlak, čas 4s"
#set title " "
#set key top right
#set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid ytics mxtics mytics lc rgb 'black' lw 0 lt 2
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0
set encoding utf8
#set xtics 0.2
set logscale y
set format y "10^{%L}"
set xrange [0:20]
set yrange [1e-8:1e-1]
set key top right
#set key at 19., 0.45
set key box maxcols 1


#set key off
#set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )


#plot bt u 1:2 w l lw 2 lc rgb "web-green" title sprintf('BT'), dbc u 1:2 w l lw 2 lc rgb "red" title sprintf('DBC'), mdbc u 1:2 w l lw 2 lc rgb "blue" title sprintf('  mDBC')
#plot full u 1:2 w l lw 2 lc rgb "red" title sprintf(' full'), lowered u 1:2 w l lw 2 lc rgb "blue" title sprintf('lowered')
plot inputDataFile u 1:2 w l lw 2 lc rgb "red" title sprintf(' E_k tot. ')

