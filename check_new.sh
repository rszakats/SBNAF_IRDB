#!/bin/bash

workdir='/var/gavo/inputs/sbnaf/data'
cd $workdir
minsize=120000000
cdate=`date`
sdate=`date +%s`
mkdir -p ${workdir}/backup
mkdir -p ${workdir}/logs

c=0
while true;do
  wget -q http://193.6.22.152/public.csv -O ${workdir}/tmp
  actualsize=$(wc -c <"tmp")
  if [ $actualsize -ge $minsize ]; then
    iscsv=`awk ' BEGIN{FS=","}!n{n=NF}n!=NF{failed=1;exit}END{print !failed}' ${workdir}/tmp`
    if [ "$iscsv" == 1 ];then
      orig=`md5sum public.csv | awk '{print $1}'`
      new=`md5sum tmp | awk '{print $1}'`
      if [ "$orig" == "$new" ];then
        echo $cdate "The new file is the same as the old. Nothing to do." >> ${workdir}/logs/log.log
      else
        echo $cdate "Replacing old file with new one." >> ${workdir}/logs/log.log
        mv public.csv backup/public.csv${sdate}
        mv tmp public.csv
        nfiles=`ls -lrt backup/ | wc -l`
        if [ "$nfiles" -ge "6" ];then
          ls -lrt backup/ | head -n 2 | tail -n 1 | awk '{print "rm backup/"$9}' | bash
        fi

        cd /var/gavo/inputs/sbnaf
        dachs imp q.rd
        gavo serve restart
      fi
      break
    else
      echo $cdate "The csv file is not valid!" >> ${workdir}/logs/error.log
      break
    fi
  else
    if [ "$c" -ge "5" ];then
      echo $cdate "Downloading the datafile failed!" >> ${workdir}/logs/error.log
      break
    fi
    c=$[$c+1]
    sleep 10
  fi
done

if [ -f "${workdir}/tmp" ]; then
  rm ${workdir}/tmp
fi

