#for i in Run02_2019_12_03_12_11_41 Run02_2019_12_05_10_31_02 Run03_2019_12_03_14_35_49 Run03_2019_12_05_11_24_42 Run04_2019_12_03_14_44_33 Run04_2019_12_05_12_48_57 Run05_2019_12_03_16_44_12 Run05_2019_12_05_13_11_05 Run06_2019_12_03_16_58_44 Run06_2019_12_05_13_35_30 Run07_2019_12_03_17_23_58 Run07_2019_12_05_13_58_46 Run08_2019_12_03_18_24_56 Run08_2019_12_05_14_37_04 Run09_2019_12_05_15_05_54 Run10_2019_12_05_15_36_22 Run11_2019_12_05_16_01_32 Run12_2019_12_05_17_00_17 Run13_2019_12_05_17_51_07 Run14_2019_12_05_18_20_41 Run15_2019_12_05_18_43_46 Run16_2019_12_05_19_08_17 Run17_2019_12_05_19_42_07;

#datapath="/media/guang/DATA1/sFGD_norm_FR_bldg17"
#datapath="/home/guang/work/sfgd_framework/data_preprocessing/LANL2020"
#datapath="/media/guang/LANL2020_1/03"
#datapath="/media/guang/LANL2020_1/comparison_old"
#datapath="/media/guang/LANL2020_1"
#datapath="/media/guang/LANL2020_1/beam/${1}/00"
datapath='/home/guang/work/testArea'

#for i in Run5_2020_11_17_11_11_09 Run6_2020_11_17_11_17_38 Run7_2020_11_17_11_20_08 Run8_2020_11_17_11_22_03;
#for i in Run0_2020_11_09_12_16_17; #Run0_2020_11_08_11_35_54 Run0_2020_11_08_11_40_13 Run0_2020_11_08_11_43_43;
#for i in Run_000_SubRun_03_2020_11_23_16_57_27;
#for i in Run0_2020_11_24_09_39_58 Run1_2020_11_24_09_40_55; 
#for i in Run_002_SubRun_00_2020_11_29_17_12_29 Run_002_SubRun_01_2020_11_29_17_17_52 Run_002_SubRun_02_2020_11_29_17_23_16 Run_002_SubRun_03_2020_11_29_17_28_39;
  for i in $1;
  do
    ./TDMunpack -f ${datapath}/MCR0_${i}.daq
    ./TDMunpack -f ${datapath}/MCR1_${i}.daq
    ./TDMunpack -f ${datapath}/MCR2_${i}.daq
    ./TDMunpack -f ${datapath}/MCR3_${i}.daq
    #./TDMunpack -f ${datapath}/MCR7_${i}.daq
    ls ${datapath}/MCR*${i}*Slot* > febs_files_list.list

    ./unpack 
    #mv temp_raw.root ${datapath}/MCR0_${i}__raw.root
    ./Calibration ${datapath}/MCR0_${i}__raw.root 55 50
    ./EventStructure ${datapath}/MCR0_${i}__calib.root LANL90
    #./EventStructure_cosmics ${datapath}/MCR0_${i}__calib.root USJ
    #  rm -rf febs_files_list.list
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
