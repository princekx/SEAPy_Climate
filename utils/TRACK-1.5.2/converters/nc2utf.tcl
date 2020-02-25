#!/users/nutis/kih/local/bin/convsh

if { $argc != 2 } {
   puts "$argv0: Usage: $argv0 <infile> <outfile>"
   exit 1
} else {
   puts "$argv"
}

set infile [lindex $argv 0]
set outfile [lindex $argv 1]

set filetype 0
set outformat utf
set fieldlist -1

readfile $filetype $infile

writefile $outformat $outfile $fieldlist
