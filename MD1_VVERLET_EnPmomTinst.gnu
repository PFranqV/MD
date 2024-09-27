#Gnuplot script: Gr▒| fic K, V, Etot vs t1' amb Verlet integrator
set autoscale
#maco
set term jpeg
set output 'MD1_VVERLET_Energs.jpeg'
set title "E'=f(t'), {/Symbol D}t'=0.0001"
set xlabel "t'"
set ylabel "E'"
set xr [0:20]
set linetype 1 lc rgb "greenyellow" lw 1 #potential#
set linetype 2 lc rgb "web-green" lw 1   #kinetic#
set linetype 3 lc rgb "web-blue" lw 1.5  #total#
#opcions
set border 3
set xtics nomirror
set ytics nomirror
set key center right
#grafic
plot "MD1termo_VVerlet.dat" every ::1 u 1:5 t"U^{LJ}'" with l lt 1, \
"MD1termo_VVerlet.dat" every ::1 u 1:4 t"K'" with l lt 2, \
"MD1termo_VVerlet.dat" every ::1 u 1:6 t"E'" with l lt 3

#Gnuplot script: Gr▒| fic K, V, Etot vs t1' amb Verlet integrator
set autoscale
set term jpeg
set output 'MD1_VVERLET_Etot.jpeg'
#maco
set title "E'=f(t'), {/Symbol D}t'=0.0001"
set xlabel "t'"
set ylabel "E'"
set xr [0:20]
set yr [18150:18220]
set linetype 1 lc rgb "web-blue" lw 1.5  #total#
#opcions
set border 3
set xtics nomirror
set ytics nomirror
set key top center
#grafic
plot "MD1termo_VVerlet.dat" every ::1 u 1:6 t"E'" with l lt 1

#Gnuplot script: Gr▒| fic K, V, Etot vs t1' amb Verlet integrator
set term jpeg
set output 'MD1_VVERLET_Pmom.jpeg'
set autoscale
set title "p'=f(t'), {/Symbol D}t'=0.0001"
set xlabel "t'"
set ylabel "p'"
set xr [0:20]
set yr [127:128]
set linetype 1 lc rgb "navy" lw 1 #pmom#
#opcions
set border 3
set xtics nomirror
set ytics nomirror
set key bottom center
#grafic
plot "MD1termo_VVerlet.dat" every ::1 u 1:2 t "p'" with l lt 1

#Gnuplot script: Gr▒| fic K, V, Etot vs t1' amb Verlet integrator
set autoscale
#maco
set term jpeg
set title "T'_{inst}=f(t'), {/Symbol D}t'=0.0001"
set output 'MD1_VVERLET_Tinst.jpeg'
set xlabel "t'"
set ylabel "T'_{inst}"
set xr [-0.5:20]
set yr [70:100]
set linetype 1 lc rgb "red" lw 1  #total#
#opcions
set border 3
set xtics nomirror
set ytics nomirror
set key top center
#grafic
plot "MD1termo_VVerlet.dat" every ::1 u 1:3 t"T'_{inst}'" with l lt 1

