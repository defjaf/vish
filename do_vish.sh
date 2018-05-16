#!/bin/sh -kuh

model=x
do_read_Pk=0
do_Gauss_vis=0
tau_r=0.4; w_r=0.5; Delta_w=0.1;
transfile="dummy";
Omega_vac=0.0;
do_XeXe=0; do_Xed=0;
Radius=10.0; z_i=20.0; dz_i=2; 

case $model in 
  x) tau_r=0.0; w_r=0.965; Delta_w=0.003;
     Omega=1; h=0.5; Omega_vac=0.0; n_s=1.0; OmegaBh2=0.0125;
     z_m=100.0; x_0=1.0; xbar=1.0; alpha=0.0;;

  knox) Omega=0.6; Omega_vac=0.2; OmegaBh2=0.025;
	z_m=26.0; x_0=1.0; tau_r=0.0; alpha=0.0
    ;;
  0) Omega=0.3; h=0.5; Omega_vac=0.0; n_s=1.0; OmegaBh2=0.0125;
     z_m=183.0; x_0=0.2; xbar=1.0; alpha=0.0; tau_r=0.0;;
  01) Omega=0.3; h=0.5; Omega_vac=0.0; n_s=1.0; OmegaBh2=0.0125;
     z_m=183.0; x_0=1.0; xbar=1.0; alpha=0.0; tau_r=0.0;;
  02) Omega=0.4; h=0.7; Omega_vac=0.0; n_s=1.0; OmegaBh2=0.0125;
     z_m=183.0; x_0=0.2; xbar=1.0; alpha=0.0; tau_r=0.0;;
  
  1)
     #Omega=1.0; h=0.5 ; n_s=-1.15; OmegaBh2=0.025
     Omega=1.0; h=0.5; n_s=1.0; OmegaBh2=0.025;
     transfile="trans_h5o10ob01_noreio.dat"
     z_m=800.0; x_0=0.1; xbar=1.0; alpha=0.0;
     ;;
  2)
     Omega=1.0; h=0.5; n_s=1.0; OmegaBh2=0.025;
     transfile="trans_h5o10ob01_reio.dat"
     z_m=800.0; x_0=42.3; xbar=1.0; alpha=-1.25;
     ;;
  3)
     Omega=1.0; h=0.5; n_s=1.0; OmegaBh2=0.0125;
     do_read_Pk=0;
     z_m=62.0; x_0=1.0; xbar=1.0; alpha=0.0;
     tau_r=0.0;
     ;;
  4)
     Omega=1.0; h=0.5; n_s=1.0; OmegaBh2=0.0125;
     do_read_Pk=0;
     z_m=183.0; x_0=0.2; xbar=1.0; alpha=0.0;
     tau_r=0.0;
     ;;
esac

#x_0=0.0; echo "setting x_0=0"

vish $Omega $h $n_s $OmegaBh2 $Omega_vac << EOF
$do_read_Pk
$transfile
$do_Gauss_vis
$z_m $x_0 $tau_r $xbar $alpha
$tau_r $w_r $Delta_w
$do_XeXe $do_Xed
$Radius $z_i $dz_i
EOF
