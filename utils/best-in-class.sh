function greppy {
  grep ^total */*/*/minimap2-endbonus.cigar-QV-table.txt | awk 'OFS="\t" {gsub("/","\t"); print $2,$3,$5,$6,$7,$8,$9,$10}'
}

echo Winner1
greppy | tableFilter.py -n 1 -s 3

echo; echo Winner2
greppy | tableFilter.py -n 1 -s 6

echo; echo ALL
greppy


echo ; echo Counts
greppy | awk '{a[$1]+=1 ; t+=1} END { for (e in a) print e"\t"a[e] ; print "total\t"t}'
