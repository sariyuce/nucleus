./tree $1 $2 >tree.txt
dot -Tpdf tree.txt -o $1"_"$2"_tree.pdf"
rm tree.txt
