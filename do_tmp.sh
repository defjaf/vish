#!/bin/sh -kuh

code=vish    #_approx

do_read_Pk=0
do_Gauss_vis=0
label="_patchy_yminmax"

w_r=0.5; Delta_w=0.1; transfile="dummy";
x_0=0.0; xbar=1.0; alpha=0.0;

do_XeXe=0; do_Xed=0;
Radius=0.001; z_i=10.0; dz_i=10.0; 

## set *one* of {tau_r, x_0, z_m}=0 to determine from the other 2

# B 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
#for model in B1 B2 B3 B4 17 18 19 29 21 22 23 24; do
#  for model in Pd1 P P1 0 Pd; do
#  for model in P1 P2 Pd1 Pd2 Pd P 0; do
  for model in P P1 P2 Pd Pd1 Pd2 P3 P4 P5 Pd3 Pd4 Pd5; do
    x_0=0.0; do_XeXe=0; do_Xed=0;
    case $model in 
    P) do_XeXe=1; do_Xed=0;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=0.001; z_i=10.0; dz_i=3.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    P1) do_XeXe=1; do_Xed=0;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=20.0; z_i=10.0; dz_i=3.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    P2) do_XeXe=1; do_Xed=0;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=1.0; z_i=10.0; dz_i=3.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;

    Pd) do_XeXe=0; do_Xed=1;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=0.001; z_i=10.0; dz_i=3.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    Pd1) do_XeXe=0; do_Xed=1;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=20.0; z_i=10.0; dz_i=3.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    Pd2) do_XeXe=0; do_Xed=1;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=1.0; z_i=10.0; dz_i=3.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;

    P3) do_XeXe=1; do_Xed=0;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=0.001; z_i=10.0; dz_i=10.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    P4) do_XeXe=1; do_Xed=0;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=20.0; z_i=10.0; dz_i=10.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    P5) do_XeXe=1; do_Xed=0;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=1.0; z_i=10.0; dz_i=10.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;

    Pd3) do_XeXe=0; do_Xed=1;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=0.001; z_i=10.0; dz_i=10.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    Pd4) do_XeXe=0; do_Xed=1;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=20.0; z_i=10.0; dz_i=10.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;
    Pd5) do_XeXe=0; do_Xed=1;
      Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125;
      Radius=1.0; z_i=10.0; dz_i=10.0; 
      z_m=5; tau_r=0.0; x_0=1.0;;

    B1) Omega=0.3; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.05;
	x_0=1.0; z_m=0; tau_r=0.1;;
    B2) Omega=0.3; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.05;
	x_0=0.2; z_m=0; tau_r=0.1;;
    B3) Omega=0.3; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.05;
	x_0=1.0; z_m=0; tau_r=0.5;;
    B4) Omega=0.3; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.05;
	x_0=0.2; z_m=0; tau_r=0.5;;

    B) Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0375; 
       z_m=56; tau_r=0.0; x_0=1.0;;
    0) Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125; 
       z_m=5; tau_r=0.0; x_0=1.0;;

    1) Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125; z_m=19; tau_r=0.1;;
    2) Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125; z_m=56; tau_r=0.1;;
    3) Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125; z_m=56; tau_r=0.5;;
    4) Omega=1; Omega_vac=0; h=0.5; n_s=1; OmegaBh2=0.0125; z_m=166; 
       tau_r=0.5;;
    
    5) Omega=1; Omega_vac=0; h=0.5; n_s=0.8; OmegaBh2=0.025; z_m=12; 
       tau_r=0.1;;
    6) Omega=1; Omega_vac=0; h=0.5; n_s=0.8; OmegaBh2=0.025; z_m=35; 
       tau_r=0.1;;
    7) Omega=1; Omega_vac=0; h=0.5; n_s=0.8; OmegaBh2=0.025; z_m=35; 
       tau_r=0.5;;
    8) Omega=1; Omega_vac=0; h=0.5; n_s=0.8; OmegaBh2=0.025; z_m=104; 
       tau_r=0.5;;
    
    9) Omega=0.4; Omega_vac=0.6; h=0.65; n_s=1; OmegaBh2=0.015; z_m=14;
       tau_r=0.1;;
    10) Omega=0.4; Omega_vac=0.6; h=0.65; n_s=1; OmegaBh2=0.015; z_m=44;
	tau_r=0.1;;
    11) Omega=0.4; Omega_vac=0.6; h=0.65; n_s=1; OmegaBh2=0.015; z_m=44;
	tau_r=0.5;;
    12) Omega=0.4; Omega_vac=0.6; h=0.65; n_s=1; OmegaBh2=0.015; z_m=130;
	tau_r=0.5;;
    
    13) Omega=1; Omega_vac=0; h=0.35; n_s=1; OmegaBh2=0.015; z_m=13;
	tau_r=0.1;;
    14) Omega=1; Omega_vac=0; h=0.35; n_s=1; OmegaBh2=0.015; z_m=39;
	tau_r=0.1;;
    15) Omega=1; Omega_vac=0; h=0.35; n_s=1; OmegaBh2=0.015; z_m=39;
	tau_r=0.5;;
    16) Omega=1; Omega_vac=0; h=0.35; n_s=1; OmegaBh2=0.015; z_m=116;
	tau_r=0.5;;
    
    17) Omega=0.4; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.0125; z_m=19;
	tau_r=0.1;;
    18) Omega=0.4; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.0125; z_m=54;
	tau_r=0.1;;
    19) Omega=0.4; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.0125; z_m=54;
	tau_r=0.5;;
    20) Omega=0.4; Omega_vac=0; h=0.7; n_s=1; OmegaBh2=0.0125; z_m=156;
	tau_r=0.5;;
    
    21) Omega=0.4; Omega_vac=0; h=0.8; n_s=1; OmegaBh2=0.0125; z_m=19;
	tau_r=0.1;;
    22) Omega=0.4; Omega_vac=0; h=0.8; n_s=1; OmegaBh2=0.0125; z_m=54;
	tau_r=0.1;;
    23) Omega=0.4; Omega_vac=0; h=0.8; n_s=1; OmegaBh2=0.0125; z_m=54;
	tau_r=0.5;;
    24) Omega=0.4; Omega_vac=0; h=0.8; n_s=1; OmegaBh2=0.0125; z_m=156;
	tau_r=0.5;;
  esac

  lis=${Omega},${Omega_vac},${h},${n_s},${OmegaBh2},${z_m},${x_0},${tau_r}
  if [ $do_XeXe -eq 1 ]; then lis=${lis},XeXe;
  fi
  if [ $do_Xed -eq 1 ]; then lis=${lis},Xed;
  fi
  if [ $do_Xed -eq 1 -o $do_XeXe -eq 1 ]; then 
    lis=${lis}_${Radius},${z_i},${dz_i};
  fi
  outdir=models${label}/${model}_${lis}
  mkdirhier $outdir
  
  echo "Doing model "${model} " = " ${lis}

  $code $Omega $h $n_s $OmegaBh2 $Omega_vac > ${outdir}/logfile  <<EOF 
$do_read_Pk
$transfile
$do_Gauss_vis
$z_m $x_0 $tau_r $xbar $alpha
$tau_r $w_r $Delta_w
$do_XeXe $do_Xed
$Radius $z_i $dz_i
EOF

  mv -f *.out ${outdir}/.

done
