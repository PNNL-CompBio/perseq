# keggftp/pathway/map.tar.gz for map.list
# keggftp/genes/fasta
# keggftp/genes/links/*

curl http://rest.kegg.jp/list/ko > $(date +"%Y%m%d")_ko_list.txt
curl http://rest.kegg.jp/list/pathway > $(date +"%Y%m%d")_pathway_list.txt