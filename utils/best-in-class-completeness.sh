function greppy {
  grep ^Completeness */*/*/miniQV.out.txt 
}

function formatty { 
  awk 'OFS="\t" {gsub("/","\t"); print $2,$3,$(NF)}' 
}

function all {
  greppy | formatty
}

function best {
  all | tableFilter.py -n 1 -s 3
}

function count {
  all | awk '{a[$1]+=1}END{for (e in a) print e"\t"a[e]}'
}

echo Best; best
echo ; echo All ; all
echo ; echo Count ; count
