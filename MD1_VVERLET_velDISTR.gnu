#Gnuplot
set autoscale
set style fill solid
set term jpeg
set output "Vxyz_VVERLET_init.jpeg"
set ylabel "Frecuencia"
set multiplot layout 1,3 title "Distribuciones por componentes de v'"

# Primer gráfico para VxA
set xlabel "v'_{x}"
set xrange [-11:11]
plot "MD1_VxyzDISTR_init.dat" u 1:2 w boxes lc rgb "blue" t "V'x"

# Segundo gráfico para Vy/
set xlabel "v'_{y}"
set xrange [-11:11]
plot "MD1_VxyzDISTR_init.dat" u 1:3 w boxes lc rgb "green" t "V'y"

# Tercer gráfico para VzA
set xlabel "v'_{z}"
set xrange [-11:11]
plot "MD1_VxyzDISTR_init.dat" u 1:4 w boxes lc rgb "red" t "V'z"

unset multiplot
# Configuración de estilo de las barras
set style fill solid border lc rgb "black"

# Vxyz_VVERLET_FINAL
set autoscale
set term jpeg
set output "Vxyz_VVERLET_fin.jpeg"
set ylabel "Frecuencia"
set multiplot layout 1,3 title "Distribuciones por componentes de v'"
set xtics font ",6"
set linestyle 1 lc rgb "brown" lw 3
set xr [-41:41]
set size ratio 1
set xlabel "Velocidad reducida, x"
plot "MD1_VxyzDISTR_fin.dat" u 1:2 w boxes lc rgb "blue" t "V'x" ,\
"MD1_VxyzDISTR_fin.dat" u 1:2 w l lc rgb "black" lw 1.5
set xr [-41:41]
set xtics font ",6"
# Segundo gráfico para Vy/
set xlabel "Velocidad reducida, y"
set size ratio 0.8
plot "MD1_VxyzDISTR_fin.dat" u 1:3 w boxes lc rgb "green" t "V'y" ,\
"MD1_VxyzDISTR_fin.dat" u 1:3 w l lc rgb "black" lw 1.5
set xr [-41:41]

# Tercer gráfico para VzA
set xlabel "Velocidad reducida, z"
set size ratio 1
plot "MD1_VxyzDISTR_fin.dat" u 1:4 w boxes lc rgb "red" t "V'z" ,\
"MD1_VxyzDISTR_fin.dat" u 1:4 w l lc rgb "black" lw 1.5
set xr [-41:41]
unset multiplot

#script for plotting distribución de las velocidades por componentes

# Vnorm_VVERLET_FINAL
set key
set term jpeg
set output "Vnorm_VVERLET_fin.jpeg"
set style fill solid border lc rgb 'light-blue'
set title "Distribución de velocidades finales"
set xlabel "|v'|"
set ylabel "Frecuencia"
plot "MD1_VnormDISTR_fin.dat" u 1:2 w boxes t "|v'|" ,\
"MD1_VnormDISTR_fin.dat" u 1:2 w l lc rgb "brown" lw 1.5


