for j in 234;
do
  for k in 0{0..9} {10..11}
  #for k in {10..11}
  do
    datapath="/media/guang/Elements/Foundfiles_LANL2020/Lostpartition2/beam/${j}/${k}"
    if [ ! -d ${datapath} ] ; then
      break
    else
      echo ${datapath}
    fi
    for file in ${datapath}/*; do
      filename="${file##*/}"
      echo ${filename}
    done

    orig=${filename}
    one=${orig#*MCR[0-9]_}
    two=${one%.*}

    printf "Result:\n"
    printf "$orig\n"
    printf "$one\n"
    printf "$two\n"

    for i in ${two};
    do
      ./TDMunpack -f ${datapath}/MCR0_${i}.daq
      ./TDMunpack -f ${datapath}/MCR1_${i}.daq
      ./TDMunpack -f ${datapath}/MCR2_${i}.daq
      ./TDMunpack -f ${datapath}/MCR3_${i}.daq
      ./TDMunpack -f ${datapath}/MCR7_${i}.daq
      ls ${datapath}/MCR*${i}*Slot* > febs_files_list.list

      ./unpack 
      #mv temp_raw.root ${datapath}/MCR0_${i}__raw.root
      ./Calibration ${datapath}/MCR0_${i}__raw.root 55 50
      ./EventStructure ${datapath}/MCR0_${i}__calib.root USJ
      #./EventStructure_cosmics ${datapath}/MCR0_${i}__calib.root USJ
    #  rm -rf febs_files_list.list
    done
  done
done    
