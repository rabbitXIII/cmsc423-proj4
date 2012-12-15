if [ $# -ne 2 ];
then
	echo "Please run as follows: $0 [SOURCE_FILE] [OUTPUT_FASTA]"
	exit
fi


SOURCE_FILE=$1
OUTPUT_FILE=$2

echo "Warning: This is not a very fast process..."

toAmos_new -s $SOURCE_FILE -b rgopal.bnk

./proj4_rgopal.pl $SOURCE_FILE HOXD2.txt rgopal_output.ovl

bank-transact -b rgopal.bnk -m rgopal_output.ovl 

tigger -b rgopal.bnk

make-consensus -B -b rgopal.bnk 

bank2fasta -b rgopal.bnk > $2

rm -rf rgopal.bnk

