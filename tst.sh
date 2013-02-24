t1=`date +%s -d "Sun Feb 24 12:20:02 MSK 2013"`
t2=`date +%s -d "20120115 21:35:48"`
diff=$(($t1 - $t2))

echo -e "  diff is\t"$(($t1 - $t2))
