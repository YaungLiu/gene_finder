for file in ncbi_dataset/ncbi_dataset/data/GCA_*/GCA_*.fna; do
  echo "Processing $file"
  python gene_finder2.py $file > $file.txt
done
