# VCF FORMAT is GT:DS:HDS:GP. This script reads the second field, dosage (DS).
# INFO has AF, MAF, AC, AN.
echo "Started 🧬🧬🧬🧬🧬🧬🧬"
export QCOPT="--filt-DP=0 --filt-GQ=0 --filt-minorVAF=0 --filt-GP=0"
export VCFDIR="/mnt/numberTwo2/Faye_GRS_2022/Genotyping/ProcessedData/Genotype_ethinicities/healthy_controlOsenbergUCSF/filteredVcfs"
export Reference="/mnt/numberSix/MirandaProjectGRS2025/Data/"
export AF=AF # allele frequency field in INFO
export DSFLD=2 # DS field in GTP
zcat $VCFDIR/healthy_controlOsenbergUCSF_GRSFayeMarch2023.vcf.gz | grep -v ^## | head -1 | cut -f 10- | tr '\t' '\n' |\
sed 's/$/\tUnknSex\tUnaff/' | sed '1 s/^/SeqID\tSex\tAff\n/' > sample_file.txt # 25

my_func() { CHR=$1; POS=$2; ID=$3; REF=$4; ALT=$5; BETA=$6; CtAF=$7 # CtAF is not used.
  if [ "$CHR" != "#CHROM" ]; then
    echo -n "Processing $ID $CHR-$POS-$REF-$ALT beta=$BETA "
   tabix -h $VCFDIR/healthy_controlOsenbergUCSF_GRSFayeMarch2023.vcf.gz $CHR:$POS-$POS | df.grep --comment=## -t 1 --if [s4=$REF][s5=$ALT] > 0.tmp
    if [ `grep -v ^# 0.tmp | wc -l` == 1 ]; then
      vQC 0.tmp $QCOPT --rv=0 --pv=0 --filt-miss-rate=1  --genome=GRCh37 --filt-no-geno=no --filt-no-var=no --spl=sample_file.txt --prevalence=0.03 >1.tmp 2>>e.tmp
      if [ `grep -v ^# 1.tmp | wc -l` == 1 ]; then
        DEFAF=`grep -v ^# 1.tmp | cut -f 8 | df.replace -f 1 --strg "value[tf[],;,$AF=]"`
        DEFDS=$(echo "2 * $DEFAF" | bc -l) # make sure all SNPs are in diploid regions
        DEFDS=`printf "%f\n" "$DEFDS"`
        df.replace 1.tmp --keep-comment -d'\t' -f 10- --strg xtibn[tf[],:,$DSFLD] | df.replace --keep-comment -d'\t' -f 10- -,/ .,$DEFDS |\
        df.replace --keep-comment -d'\t' -f 10- --math "tf[]*$BETA-$DEFDS*$BETA" > 2.tmp
        if [ -s DSxBETA.txt ]; then grep -v ^# 2.tmp >> DSxBETA.txt
        else cat 2.tmp > DSxBETA.txt
        fi
        echo done
      else
        REASON=`grep ':          1 ' e.tmp | sed 's/^.* //' | tail -1`
        echo "failed: filtered out due to $REASON"
      fi
    else
      echo "failed: not in the imputed VCF"
    fi
  fi ; }
export -f my_func
rm -f *.tmp DSxBETA.txt*
parallel -j 1 -a $Reference/FinalSNP_TsoiSupp4_AB_04292022.vcf --colsep '\t' my_func
# copy screen output to DS_prs.all.log
gzip -f DSxBETA.txt

for SNPSET in all noHLA HLA; do
  if   [ "$SNPSET" == all ]; then PREFIX=$Reference/FinalSNP_TsoiSupp4_AB_04292022
  elif [ "$SNPSET" == noHLA ]; then PREFIX=$Reference/FinalSNP_Tsoi_PsO_NoMHC_Alleles_2025
  elif [ "$SNPSET" == HLA ]; then PREFIX=$Reference/FinalSNP_Tsoi_PsO_MHC_Alleles_2025
  fi
  df.grep DSxBETA.txt.gz --comment=## -t 1 --if s1,2,4,5km1,2,4,5:$PREFIX.vcf |\
  cut -f 3,10- | df.transpose | df.insert --titles PRS --math nm_stat[SUM,2-] | df.cut -f 1,-1 |\
  sed '1 s/ID/StudyID/' |\
  df.sort -t 1 -c n1a > DS_prs.$SNPSET.txt
  df.replace DS_prs.$SNPSET.txt -t 1 -f 2 --math exp[tf[]] | sed '1 s/PRS/PRS_OR/' > DS_por.$SNPSET.txt
done

rm -f *.tmp *.gz

echo "✅ Polygenic Risk Calculation has done "
echo "\n ⚡ File Saved: ${pwd}"

