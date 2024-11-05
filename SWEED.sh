#用SWEED计算CLR
#/data01/wangyf/software/sweed/SweeD -h

#拆分成50条染色体的vcf
for i in {572..621}; do
vcftools --gzvcf ../test.vcf.gz --recode --recode-INFO-all --stdout --chr NC_056${i}.1 > NC_056${i}.1.vcf;
done

#将每个chr拆成100kb的窗口，也就是length/100000，计算grid参数大小，保存到文件chr_gird.txt中

#运行SWEED
for line in `cat chr_grid.txt`
do
  chr=${line%%,*}
  grid=${line##*,}
  /data01/wangyf/software/sweed/SweeD -name chr${chr}.100kb \
        -input ${chr}.vcf \
        -grid ${grid} \
        -minsnps 40 \
        -maf 0.05 \
        -missing 0.1
sed -i '1,3d' SweeD_Report.chr${chr}.100kb #删除每个report文件的前3行,"Position Likelihood Alpha StartPos EndPos"
done

#添加染色体号
for line in `cat chr_grid.txt`
do
  chr=${line%%,*}
  awk -v chr="$chr" '{print chr,$4,$5,$3}' SweeD_Report.chr${chr}.100kb > SweeD_Report.chr${chr}.100kb-2
  sed -i "s@ @\t@g" SweeD_Report.chr${chr}.100kb-2
done

#合并
cat *-2 >total_chr.txt
awk '{sub(/\.0000/, "", $2); sub(/\.0000/, "", $3); print}' total_chr.txt  > total_chr.txt-1
sed -i "s@ @\t@g" total_chr.txt-1

#对第3列进行排序，提出前5%
awk '$3 ~ /e\+02/' total_chr.txt |sort -nr -k3 > sort.total_chr.txt
awk '$3 ~ /e\+01/' total_chr.txt |sort -nr -k3 >> sort.total_chr.txt
awk '$3 ~ /e\+00/' total_chr.txt |sort -nr -k3 |head -494 >> sort.total_chr.txt

#前1%,提出第1、5、6列形成.bed文件
head -152 sort.total_chr.txt |awk '{print $1 "\t" $5 "\t" $6}' > sort1%.total_chr.txt

#R可视化
conda activate r4.3

chr1 <- read.table(file="SweeD_Report.NC_056572.1.100kb",sep = "\t", header = T)
plot(chr1[,1]/10000,chr1[,2])
