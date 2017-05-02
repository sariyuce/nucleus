./ktruss_dtruss3 $1 1 0 >$1"_1_allsubs"
./ktruss_dtruss3 $1 2 0 >$1"_2_allsubs"
./ktruss_dtruss3 $1 3 0 >$1"_3_allsubs"

cat $1"_1_allsubs"|grep density|awk '-F[,:()]' '{print $1" "$3" "$6" "$7" "$12" "$14}' >$1"_1_allsubs_res"
cat $1"_2_allsubs"|grep density|awk '-F[,:()]' '{print $1" "$3" "$6" "$7" "$12" "$14}' >$1"_2_allsubs_res"
cat $1"_3_allsubs"|grep density|awk '-F[,:()]' '{print $1" "$3" "$6" "$7" "$12" "$14}' >$1"_3_allsubs_res"
