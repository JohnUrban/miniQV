function greppy {
  grep ^Completeness */*/*/miniQV.out.txt 
}

function formatty { 
  awk 'OFS="\t" {gsub("/","\t"); print $2,$3,$(NF)}' 
}

function all {
  greppy | formatty
}

all
