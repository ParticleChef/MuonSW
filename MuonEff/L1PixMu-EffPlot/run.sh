rm -rf eff_plot_v3_C.so
rm -rf eff_plot_v3_C.d
rm -rf result.log

root -l -b  < x_file.C  &> result.log &

rm -rf eff_plot_v3_C.so
rm -rf eff_plot_v3_C.d
