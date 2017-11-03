EXT=fasta
for i in *; do
  if [ "${i}" != "${i%.${EXT}}" ];then
    echo "processing $i"
    python2 dynamine.py $i
  fi
done
