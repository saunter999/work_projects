reset
set   autoscale                        # scale axes automatically
set term postscript eps color blacktext "Times" 26
#set logscale x 
#set logscale y 
set output 'phase.eps'
set ylabel 'M'
#set log y
#set ylabel rotate right
set xlabel 'fai'
set xtics font "Times, 12"
set ytics font "Times, 12"
set title 'Haldane model with interacion'
set key box linestyle 1 lc rgb "black"
set key Left left samplen 0.1 reverse  
set key spacing 0.6 font "Times, 16" 
#set key width -3.0
set pointsize 0.2
set style line 1 lw 1.0 lc rgb "forest-green" lt 1
set style line 2 lw 1.0 lc rgb "tan1" lt 1

plot	'v1_plot' using 1:2 notitle with linespoints pt 7 lt 3 lc rgb "black",\
        sin(x)*3*1.732,\
        -sin(x)*3*1.732
                    

