#!/bin/sh -kuh

label="_bdm"
do_read_Pk=1
do_Gauss_vis=0
tau_r=0.4; w_r=0.5; Delta_w=0.1;
transfile="dummy";
Omega_vac=0.0;
do_XeXe=0; do_Xed=0;
Radius=10.0; z_i=20.0; dz_i=2; 
tau_r=0.0;

for model in 4 5 4r 5r; do
case $model in 
  1) transfile="bdm/trans_h7o03ob01_noreio.dat";
     Omega=0.3; Omega_vac=0; h=0.7; n_s=2.85; OmegaBh2=0.049;
     x_0=0.1; xbar=0.1; alpha=0.0; z_m=800.0
    ;;
  2) transfile="bdm/trans_h5o10ob01_reio.dat"
     Omega=1; Omega_vac=0; h=0.5; n_s=2.85; OmegaBh2=0.025;
     x_0=42.3; xbar=1.0; alpha=-1.24; z_m=800.0
    ;;
  3) transfile="bdm/trans_h5o10ob01_noreio.dat"
     Omega=1; Omega_vac=0; h=0.5; n_s=2.85; OmegaBh2=0.025;
     x_0=0.1; xbar=0.1; alpha=0.0; z_m=800.0
    ;;
  4r) transfile="bdm/transh7o02ob05.dat";
      Omega=0.2; Omega_vac=0; h=0.7; n_s=2.25.; OmegaBh2=0.0245;
     x_0=42.3; xbar=1.0; alpha=-1.24; z_m=800.0
    ;;
  5r) transfile="bdm/transh7o04ob10.dat";
     Omega=0.4; Omega_vac=0; h=0.7; n_s=2.05; OmegaBh2=0.049;
     x_0=42.3; xbar=1.0; alpha=-1.24; z_m=800.0
    ;;  
  4) transfile="bdm/transh7o02ob05.dat";
     Omega=0.2; Omega_vac=0; h=0.7; n_s=2.25; OmegaBh2=0.0245;
     x_0=0.1; xbar=0.1; alpha=0.0; z_m=800.0
    ;;
  5) transfile="bdm/transh7o04ob10.dat";
     Omega=0.4; Omega_vac=0; h=0.7; n_s=2.05; OmegaBh2=0.045;
     x_0=0.1; xbar=0.1; alpha=0.0; z_m=800.0
     ;;  
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

vish $Omega $h $n_s $OmegaBh2 $Omega_vac  > ${outdir}/logfile << EOF
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
