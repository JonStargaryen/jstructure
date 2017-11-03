EXT=pdb
for i in *; do
  if [ "${i}" != "${i%.${EXT}}" ];then
    pdb=${i%.*}
    echo "processing $pdb"
    echo "python2 Kink_Finder.py -f /home/bittrich/Downloads/KF_err_lin/$pdb.pdb -o results_$pdb/ -d"
    python2 Kink_Finder.py -f /home/bittrich/Downloads/KF_err_lin/$pdb.pdb -o results_$pdb/ -d
  fi
done
