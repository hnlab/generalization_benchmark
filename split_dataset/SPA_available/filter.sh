file="/home/zhuhui/dataset/general_refine/pdbbind_2020.code"
for pdb in `cat ${file}`
do 
    protein="/home/zhuhui/dataset/general_refine/${pdb}/${pdb}_protein.pdb"
    echo ${pdb} >> ./out
    cat ${protein} |grep HETATM |cut -d " " -f 3-4 |uniq |xargs >> ./out
done