#!/bin/bash
# =======================================================================
# files nad params
# =======================================================================

# set vars and files
yearStart=1870
yearEnd=2005
inDir="/data6/okazaki/MIROC5.0/out"
expName="historical_MIROC5_r2i1p1"
#expName="historical_MIROC5_r1i1p1"
#expName="historical_hadisst_re"
hres="t42"
nlev=18
tmpDir="../data/work"
outDir="../data/${expName}"
outFixDir="../data/fix/${hres}"

# input files (org)
gtDir="${inDir}/${expName}"
vars_gt2pb="prcp prcp01 evap evap01 T2 rh2 Ps qt qt01"
xxFileSST="/data6/okazaki/CMIP5/historical/gt/tos_Omon_MIROC5_historical_r2i1p1_t42_185001-201212.gt"
#xxFileSST="/data6/okazaki/CMIP5/historical/gt/tos_Omon_MIROC5_historical_r1i1p1_t42_185001-201212.gt"
sstYearStart=1850
sstYearEnd=2012
#xxFileSST="/data11/okazaki/HadISST/data/nc/202104/HadISST_sst.nc"
#sstYearStart=1870
#sstYearEnd=2020
ncFileSSWIso="/data6/okazaki/LeGrande/calculated_d18O_v1_1.nc"
ncFileSDWIso="/data6/okazaki/LeGrande/calculated_d18O_v1_1.nc"
ctlFileNPP="/home/okazaki/Ref/NPP/grid_npp.ctl"
refFile="/home/okazaki/Ref/CTL/${hres}.ctl"

# check
if [[ ${expName} == *MIROC* && ${xxFileSST} != *MIROC* ]] ; then
    echo 'u better off check sst file.'
fi
if [[ ${expName} == *MIROC* && ${sstYearStart} != 1850 ]] ; then
    echo 'u better off check sstYearStart & sstYearEnd'
fi

# fortran compile flags
fflag="-O -assume byterecl -convert big_endian -mcmodel=large -shared-intel"

# =======================================================================
# proxy system models
# =======================================================================

# run = 1, no = 0
#     icecore coral treering speleothem
#-----------------------------------------------
runPSM=(    1     1        1          1)
#runPSM=(    0     0        0          1)
## ice core
runIcecore=${runPSM[0]}
psmIcecore="psm_icecore.f90"
#psmIcecore="psm_icecore_dee2015.f90"
psmIcecoreX="psmIcecore.x"
## coral
runCoral=${runPSM[1]}
psmCoral="psm_coral_liu2013.f90"
psmCoralX="psmCcoral.x"
## tree-ring cellulose
runTreering=${runPSM[2]}
psmTreering="psm_treering-cellulose_roden2000.f90"
psmTreeringX="psmTreering.x"
## speleothem
runSpeleothem=${runPSM[3]}
#psmSpeleothem="psm_speleothem_bakerbradley2010.f90"
psmSpeleothem="psm_speleothem_dee2015.f90"
psmSpeleothemX="psmSpeleothem.x"


# =======================================================================
# internal & output files (automatically defined in default)
# =======================================================================

# input files (plain binary)
pbFilePrcp="${tmpDir}/prcp.bin"
pbFilePrcpIso="${tmpDir}/prcp01.bin"
pbFileEvap="${tmpDir}/evap.bin"
pbFileEvapIso="${tmpDir}/evap01.bin"
pbFileQv="${tmpDir}/qt.bin"
pbFileQvIso="${tmpDir}/qt01.bin"
pbFileT2="${tmpDir}/T2.bin"
pbFileRH2="${tmpDir}/rh2.bin"
pbFilePs="${tmpDir}/Ps.bin"
suf=$(echo $xxFileSST | rev | cut -d "." -f 1 | rev)
pbFileSST="${outFixDir}/$(basename $xxFileSST ."${suf}")_${yearStart}-${yearEnd}.bin"
pbFileSSWIso="${outFixDir}/$(basename $ncFileSSWIso .nc)_z33.bin"
pbFileSDWIso="${outFixDir}/$(basename $ncFileSDWIso .nc)_z29.bin"
pbFileNPP="${outFixDir}/$(basename $ctlFileNPP .ctl).bin"


# output file (plain binary)
pbFileIcecore="${outDir}/$(basename $psmIcecore .f90 | cut -c 5-)_${yearStart}-${yearEnd}.bin"
pbFileCoral="${outDir}/$(basename $psmCoral .f90 | cut -c 5-)_${yearStart}-${yearEnd}.bin"
pbFileTreering="${outDir}/$(basename $psmTreering .f90 | cut -c 5-)_${yearStart}-${yearEnd}.bin"
pbFileSpeleothem="${outDir}/$(basename $psmSpeleothem .f90 | cut -c 5-)_${yearStart}-${yearEnd}.bin"


# model resolution
case ${hres} in
    "t42") nx=128 ; ny=64  ; ijmax=$((nx*ny)) ;;
    "t63") nx=192 ; ny=94  ; ijmax=$((nx*ny)) ;;
    "t85") nx=256 ; ny=128 ; ijmax=$((nx*ny)) ;;
    *) echo "horizontal resolution $hres is not supported. --Aborted." ; exit ;;
esac
nt=$(((yearEnd - yearStart + 1) * 12))

# create directory
if [ ! -d ${tmpDir} ] ; then mkdir -p ${tmpDir} ; fi
if [ ! -d ${outDir} ] ; then mkdir -p ${outDir} ; fi
if [ ! -d ${outFixDir} ] ; then mkdir -p ${outFixDir} ; fi


# =======================================================================
# prepare input files for PSM
# =======================================================================

## convert GTOOL to plain binary
for var in ${vars_gt2pb} ; do
    pbFile="${tmpDir}/${var}.bin"
    if [ ! -s ${pbFile} ] ; then
        # check file
        if [ ! -s ${gtDir}/y${yearStart}/${var} ] ; then echo "File not found : ${gtDir}/y${yearStart}/${var}" ; exit ; fi
        # concatenate GTOOL files
        ./cat.sh ${gtDir} ${var} ${yearStart} ${yearEnd} ${tmpDir}/${var}
        # convert GTOOL to plain binary
        if [ ${var} = "qt" -o ${var} = "qt01" ] ; then nz=18 ; else nz=1 ; fi
        offset=0
        fact=1
        gradsgt -blc "gt2pb.gs ${tmpDir}/${var} ${pbFile} ${var} ${nx} ${ny} 1 ${nz} 1 ${nt} ${offset} ${fact}"
        # rm gtool data
        rm ${tmpDir}/${var}
    fi
done
## convert SST file (nc/gt) to binary
if [ ! -s ${pbFileSST} ] ; then
    if [ ${yearStart} -lt ${sstYearStart} -o ${yearStart} -gt ${sstYearEnd} -o ${yearEnd} -lt ${sstYearStart} -o ${yearEnd} -gt ${sstYearEnd} ] ; then
        echo "targeting years exceed periof of HadISST. --Aborted." ; exit
    fi
    tstr=$(((yearStart - sstYearStart) * 12 + 1))
    tend=$(((yearEnd - sstYearStart + 1) * 12))
    suf=$(echo $xxFileSST | rev | cut -d "." -f 1 | rev)
    case ${suf} in
    "gt")
        # MIROC
        offset=-273.15
        fact=1
        gradsgt -blc "gt2pb.gs ${xxFileSST} ${pbFileSST} tos_omon_miroc5 ${nx} ${ny} 1 1 ${tstr} ${tend} ${offset} ${fact}"
        ;;
    "nc")
        # HadISST
        grads2.1 -blc "nc2pb.gs ${xxFileSST} ${pbFileSST} ${refFile} sst ref ${nx} ${ny} 1 1 ${tstr} ${tend}"
        ;;
    esac
fi
## convert NetCDF to plain binary
#### Sea water isotope (LeGrande and Schmidt, 2006)
if [ ! -s ${pbFileSSWIso} -o ! -s ${pbFileSDWIso} ] ; then
    grads2.1 -blc "nc2pb.gs ${ncFileSSWIso} ${pbFileSSWIso} ${refFile} d18o ref ${nx} ${ny} 33 33 1 1"
    grads2.1 -blc "nc2pb.gs ${ncFileSDWIso} ${pbFileSDWIso} ${refFile} d18o ref ${nx} ${ny} 29 29 1 1"
fi
#### Net Primary Production  
if [ ! -s ${pbFileNPP} ] ; then
    grads2.1 -blc "grd2pb.gs ${ctlFileNPP} ${pbFileNPP} ${refFile} npp ref ${nx} ${ny} 1 1 1 12"
fi

# =======================================================================
# compile and run PSM
# =======================================================================

## ice core
if [ ${runIcecore} -eq 1 ] ; then
    cat ${psmIcecore} |\
    sed -e "s/<yearStart>/${yearStart}/g" \
        -e "s/<yearEnd>/${yearEnd}/g" \
        -e "s/<ijmax>/${ijmax}/g" \
    > tmp.f90
    ifort ${fflag} tmp.f90 -o ${psmIcecoreX}
    case $(basename $psmIcecore .f90) in
    "psm_icecore") 
        ./${psmIcecoreX} ${pbFilePrcp} ${pbFilePrcpIso} \
                    ${pbFileIcecore}
        ;;
    "psm_icecore_dee2015") 
        ./${psmIcecoreX} ${pbFilePrcp} ${pbFilePrcpIso} ${pbFileT2} ${pbFilePs} \
                    ${pbFileIcecore}
        ;;
    esac
fi
## coral
if [ ${runCoral} -eq 1 ] ; then
    cat ${psmCoral} |\
    sed -e "s/<yearStart>/${yearStart}/g" \
        -e "s/<yearEnd>/${yearEnd}/g" \
        -e "s/<ijmax>/${ijmax}/g" \
    > tmp.f90
    ifort ${fflag} tmp.f90 -o ${psmCoralX}
    case $(basename $psmCoral .f90) in
    "psm_coral_liu2013") 
        ./${psmCoralX} ${pbFilePrcp} ${pbFilePrcpIso} ${pbFileEvap} ${pbFileEvapIso} \
                   ${pbFileSST}  ${pbFileSSWIso}  ${pbFileSDWIso} \
                   ${pbFileCoral}
        ;;
    esac
fi
## tree-ring cellulose
if [ ${runTreering} -eq 1 ] ; then
    cat ${psmTreering} |\
    sed -e "s/<yearStart>/${yearStart}/g" \
        -e "s/<yearEnd>/${yearEnd}/g" \
        -e "s/<ijmax>/${ijmax}/g" \
        -e "s/<nz>/${nlev}/g" \
    > tmp.f90
    ifort ${fflag} tmp.f90 -o ${psmTreeringX}
    case $(basename $psmTreering .f90) in
    "psm_treering-cellulose_roden2000")
        ./${psmTreeringX} ${pbFilePrcp} ${pbFilePrcpIso} ${pbFileQv} ${pbFileQvIso} \
                    ${pbFileT2} ${pbFileRH2} ${pbFilePs} ${pbFileNPP} \
                    ${pbFileTreering}
        ;;
    esac
fi
## speleothem
if [ ${runSpeleothem} -eq 1 ] ; then
    cat ${psmSpeleothem} |\
    sed -e "s/<yearStart>/${yearStart}/g" \
        -e "s/<yearEnd>/${yearEnd}/g" \
        -e "s/<ijmax>/${ijmax}/g" \
    > tmp.f90
    ifort ${fflag} tmp.f90 -o ${psmSpeleothemX}
    case $(basename $psmSpeleothem .f90) in
    "psm_speleothem_dee2015")
        ./${psmSpeleothemX} ${pbFilePrcp} ${pbFilePrcpIso} ${pbFileEvap} ${pbFileEvapIso} \
                    ${pbFileT2} \
                    ${pbFileSpeleothem} \
                    ${option}
        ;;
    "psm_speleothem_bakerbradley2010")
        ./${psmSpeleothemX} ${pbFilePrcp} ${pbFilePrcpIso} ${pbFileEvap} ${pbFileEvapIso} \
                    ${pbFileT2} \
                    ${pbFileSpeleothem} \
                    ${option}
        ;;
    esac
fi

# =======================================================================
# clean up
# =======================================================================

rm ${tmpDir}/*.bin
