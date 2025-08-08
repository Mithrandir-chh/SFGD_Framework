#for i in Run02_2019_12_03_12_11_41 Run02_2019_12_05_10_31_02 Run03_2019_12_03_14_35_49 Run03_2019_12_05_11_24_42 Run04_2019_12_03_14_44_33 Run04_2019_12_05_12_48_57 Run05_2019_12_03_16_44_12 Run05_2019_12_05_13_11_05 Run06_2019_12_03_16_58_44 Run06_2019_12_05_13_35_30 Run07_2019_12_03_17_23_58 Run07_2019_12_05_13_58_46 Run08_2019_12_03_18_24_56 Run08_2019_12_05_14_37_04 Run09_2019_12_05_15_05_54 Run10_2019_12_05_15_36_22 Run11_2019_12_05_16_01_32 Run12_2019_12_05_17_00_17 Run13_2019_12_05_17_51_07 Run14_2019_12_05_18_20_41 Run15_2019_12_05_18_43_46 Run16_2019_12_05_19_08_17 Run17_2019_12_05_19_42_07;

#datapath="/media/guang/DATA1/sFGD_norm_FR_bldg17"
#datapath="/home/guang/work/sfgd_framework/data_preprocessing/LANL2020"
#datapath="/media/guang/LANL2020_1/03"
#datapath="/media/guang/LANL2020_1/comparison_old"
#datapath="/media/guang/LANL2020_1"

# finished 270..273; 
# finished 243..246;
for j in 234;
do
  for k in 0{0..9} {10..11}
  #for k in {10..11}
  do
    datapath="/media/guang/Elements/Foundfiles_LANL2020/Lostpartition2/beam/${j}/${k}"
    #datapath="/media/guang/LANL2020_1/beam/${j}/${k}"
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
#for i in Run_009_SubRun_00_2020_12_01_15_59_58;
#do
#  ./TDMunpack -f ${datapath}/MCR7_${i}.daq
#
#  ls ${datapath}/MCR*${i}*Slot* > febs_files_list.list
#
#  ./unpack
#  rm -rf febs_files_list.list
#done
