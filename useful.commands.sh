#!/bin/bash

##########################################################################################
## Bash

# To use local/path version of python/R/etc at top of script
#!/usr/bin/env python
#!/usr/bin/env Rscript

filename=$(basename "$fullfile")
extension="${filename##*.}"
filename_no_ext="${filename%.*}"

# Pass a tab on the command line
python test.py $'\t'

# Make a directory; -p makes parent dirs as needed
if [ ! -d "path/dir" ]; then
    mkdir -p path/dir
fi

# translate characters uppercase lowercase
echo -e "A\nB\nC" | tr '[:upper:]' '[:lower:]'
echo -e "A\nB\nC" | perl -pe 'chomp if eof' | tr '\n' ','
echo -ne "test\tX" > test.txt # no newline at end

# iterate through array in bash
arr=( one two three )
for i in "${arr[@]}"
do
   echo "$i" # or do whatever with individual element of the array
done

for i in `seq 1 10`;
do
    echo $i
done

for i in {8..20..2};
do
    echo $i
done

# For each line in a file
while read p; do
  echo $p
done < file.txt

# Splitting into variables
while IFS=" " read -r value1 value2 remainder; do
    ...
done < "input.txt"

# print the header (the first line of input)
# and then run the specified command on the body (the rest of the input)
# use it in a pipeline, e.g. ps | body grep somepattern
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
# Example:
ps -o pid,comm | body sort -k2


# count unique values - note that sort is critical for uniq to work
cat file.txt | cut -f 2 | sort | uniq | wc -l

tar -zcvf newarchive.tar.gz dirToArchive
tar -ztf archive.tar.gz [file to list]

ln -s targetfile linkname

# arithmetic in bash: runs bjobs and counts lines, and does subtraction using $((A-B))
alias bjc='echo $((`bjobs | wc -l` - 1))'
z=`expr $z + 3` # The 'expr' command performs the expansion.
z=$((z+3)) # The use of backticks (backquotes) in arithmetic expansion has been superseded by double parentheses -- ((...)) and $((...)) -- and also by the very convenient let construction.

# Get file lines starting at number 123
tail -n +123 file.txt

# <(cmd) is bash process substitution, where the output of the command can be treated
# as a file descriptor input for standard unix commands
paste <(cut -f 1 --complement file1.txt) <(cut -f 5 file2.txt)

sed '1d' file.txt | awk 'BEGIN {OFS="\t"}{if ($1 ~ /NA/) $1=0.5; print $1,$2,"chr9",$2"_"$3}'

# Sum fields in a file
cat file.txt | awk '{ if (NR>1){ sum+=$2 }} END {print sum}'

awk 'BEGIN {OFS="\t"} function abs(x){return ((x < 0.0) ? -x : x)} {print abs($2-$1)}'

# Pass a variable to awk
variable="line one\nline two"
awk -v var="$variable" 'BEGIN {print var}'

# Get non-header rows of VCF
zcat file.vcf.gz | grep -v "#"

# Rename a set of files
for f in DIR/VORX.*; do
    NAME=`echo $f | perl -ne '@fparts=split(/\./); print join(".", @fparts[1..$#fparts]);'`
    mv VORX.$NAME ZUTA.$NAME
done

# Get VCF samples
zcat file.vcf.gz | head -n 1000 | grep "#CHROM" | cut -f 10- | tr '\t' '\n' > vcfheader.samples.txt

# Get autosome VCF coords
zcat file.vcf.gz | grep -v "#" | perl -ne '@l=split();($l[0] =~ /^[\d]+$/) and print;' | cut -f 1,2,3 > file.snp_coords.txt

# Find files and delete them
find . -name '*.txt' -delete
# Find files and get their total size
find . -name "*.bg" -print0 | du --files0-from=- -hc | tail -n1

# Find files and grep them for 'Assigned'
find . -name '*.sh' -print0 | xargs -r0 grep -H 'Assigned'
# Find in files matching multiple patterns
find . \( -iname \*.sh -o -iname \*.py -o -iname \*.R -o -iname \*.Rmd \) -print0 | xargs -r0 grep -H 'tofind'

# Find files and send them to tail
find . -name "*.ASEcounts" | xargs -I {} sh -c "echo {}; tail -n2 {}"
# Find files and send them to gzip
find . -name '*.txt' -print0 | xargs -r0 gzip
# Find files changed in the last day
find . -name '*.txt' -mtime -1
# Find files changed in the last hour
find . -name '*.txt' -mmin -60

# Get file sizes sorted by size
du . -h -- * | sort -rh | less


bcftools view --samples-file list.txt vcf.gz | bcftools filter -O z -i 'MAF[0] >= 0.05' > maf0.05.vcf.gz

# Efficient way to extract specific SNPs from a VCF
grep -wFf rsid.list <(gunzip -c vcf.gz)

# Regular expression matching a line that does NOT contain a string
^((?!mystring).)*$\n

# Compare gzipped files to see if they are identical
cmp -i 8 file1.gz file2.gz # fastest, see http://unix.stackexchange.com/questions/64200/how-can-i-check-if-two-gzipped-files-are-equal
zdiff file1.gz file2.gz # uncompresses to compare actual contents

# Check if file exists (and is a regular file, not a device file)
if [ -f $FILE ]; then
fi

# -a archive (recurse dirs, copy symlinks as links, preserve permissions, file times, groups)
# -u update (do not replace if file newer at dest)
rsync -a -u src dest

# Change group ownership for files recursively
chown -R js29:newgroup * .

# Set file access control list to allow a user not in your group to access files
# May need to set access for directories *above* the one you want to give access to
setfacl -R -m u:nk5:rx /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/
setfacl -R -m g:team170:rx /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/

##########################################################################################
## Perl

# -ne runs per line; -a autosplits $_ into @F; -s allows variable passing with -- after script
cat file.txt | perl -sane 'chomp; print join("\t", @F, $locus)."\n"' -- -locus=BMI_2

# extract field from regex
if ($str =~ /(regex)/) {
  my $match = $1;
} else {die "match failed";}

my ($a, $b) = split /:/, $str;

##########################################################################################
## Python

http://mediawiki.internal.sanger.ac.uk/index.php/Python_virtualenv
# Set up a virtual environment for python (e.g. in a directory called python 2.7)
virtualenv $HOME/python2.7
# Use that copy of python
source $HOME/python2.7/bin/activate
# Stop using this python virtualenv
deactivate


##########################################################################################
## R 

# Get the size of all objects in memory
sort( sapply(ls(),function(x){object.size(get(x))}))

# Limit numbers to certain precision in output
rpkm[,1:ncol] <- sprintf("%.4f", rpkm[,1:ncol])
rpkm[as.numeric(rpkm) == 0] <- "0"

# ggplot rotated axis text
ggplot()... + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggplot remove legend
scale_colour_discrete(guide=F)
scale_color_manual(guide=F, values=colorScale)
scale_shape_discrete(guide=F)
# legend top right
d <- ggplot(mtcars, aes(x=wt, y=mpg, colour=cyl)) + geom_point(aes(colour=cyl)) + labs(title = "Legend is top right") +
theme(legend.justification = c(1, 1), legend.position = c(1, 1))


# Install old version of a package
require(devtools)
install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")

ViewDups = function(df, col) {
  df = as.data.frame(df)
  dupVals = df[duplicated(df[,col]), col]
  dups = df[df[,col] %in% dupVals,]
  View(dups[order(dups[,col]),])
}


##########################################################################################
## git 
# https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging

git add <filename> # stage files for commit
git status # see staged files
git commit -m <message>
git push origin master
git diff origin/master [local-path]

git branch hotfix # make a branch
git checkout hotfix # switch HEAD to point to latest branch
git checkout -b hotfix # Make hotfix branch and switch to it
# Do commits to hotfix

git checkout master
git merge hotfix # Master and hotfix now point to the same place
git mergetool # Useful if there are conflicts to merge
git branch -d hotfix # Delete branch no longer needed (i.e. after it's merged into master)



git log --oneline --decorate # see commits and branches

# github.com/jeremy37 token:
# 9679f5feaff127df93c55ac0b40638d4db915218



##########################################################################################
## Farm

# Check your priority on a queue
bqueues -r queuename

# What's occupying the yesterday queue
bjobs -u all -q yesterday


# Example of getting CPU runtime from many farm output files
grep CPU FarmOut/*.txt > runRasqual.CPUtimes.txt
cat runRasqual.CPUtimes.txt | perl -ne '@l=split(/\.|_|(\s)+/);print join("\t", @l)' | cut -f 3,17 | sort -nk1,1 > runRasqual.CPUtimes.sorted.txt

kinit
imeta qu -z seq -d study = 'Mapping regulatory variation in sensory neurons using IPS lines from the HIPSCI project' and target = 1 and manual_qc = 1 > irods.study.meta.txt

lfs quota /lustre/scratch109

#Use lfs quota for yourself and groups you are in:
lfs quota -h (-g GROUP | -u USER) FILESYSTEM

lfs quota -g otcoregen /lustre/scratch115


#HGI's LustreTree webapp also shows this (and more)
#information for all directories on Lustre:
  https://hgi.dev.sanger.ac.uk/lustretree/


# Use bash process substitution with bsub
#1. Workaround, e.g. using a pipe into grep:
bsub -o getSNPs.txt -J getSNPs "gunzip -c my.vcf.gz | grep -wFf SNPids.txt > my.snps.vcf"

#2. using bash -c, e.g.:
bsub -o getSNPs.txt -J getSNPs "/bin/bash -c 'grep -wFf SNPids.txt <(gunzip -c my.vcf.gz) > my.snps.vcf'"

# Waiting for a job to finish
submitJobs.py --MEM 200 -j jobname -o farmOut --blocking -c "sleep 5"
bsub -K -o out.1 sleep 10 &
bsub -K -o out.2 sleep 5 &
wait


##########################################################################################
## tmux
https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/


tmux		# start session
tmux ls		# list existing sessions
tmux attach -t 0	#

# Type Ctrl-b and release, then next key
# Ctrl-b %			split pane
# Ctrl-b arrow		switch panes
# Ctrl-d or exit	close pane
# Ctrl-b c			create window
# Ctrl-b p/n		previous/next window
# Ctrl-b <number>	go to window number
# Ctrl-b d			detach


##########################################################################################
## Google cloud

gcloud auth application-default login

# specify gsutil project, e.g. for copy
gsutil -u "bill-this-project" cp src dest

# docs
https://cloud.google.com/storage/docs/gsutil/commands/cp

# Useful options
# -n  no-clobber (don't replace existing items at destination)
# -r  recursive
# -c  continue if one file has an error
# -z  compress for upload (but actual files are left uncompressed)
# -P  preserve attributes (e.g. mod time, owner, group, etc)

# copy without over-writing (-n)
gsutil -m cp -r -n -L manifest_log_file dir gs://my-bucket

# When using -m, I got this error:
# Reauthentication challenge could not be answered because you are not in an interactive session.
# I couldn't figure it out, until I removed the -m option... then I was able to
# enter my password, and then all subsequent commands with -m worked.

# Get file sizes
# -c to include total size, -h for human readable, 
gsutil -m du -ch gs://bucket


# ssh to a cloud VM instance
ssh -i /Users/jeremys/.ssh/gcloud-ssh js29@104.155.63.97



