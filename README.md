# MERD: Meiotic Recombination Detection Pipeline with PacBio HiFi sequencing reads 
## Introduction

MERD is a pipeline for meiotic recombination detection with long read sequencing data. It detects meiotic events with 
trio samples. It now supports HiFi sequencing data only. We obtained HIGH-CONFIDENCE Single Nucleotide Variants (SNVs)
with HiFi data. The SNVs from the family are phased with both their genetic maps and HiFi reads. Then the recombination
signals are extracted according to the shift error two groups of genetic maps. Finally, We filter and cluster the extracted 
signals to obtain the meiotic Recombination events.

## Quick start 
### 1. Prepare and install required software and package 
* pysam
* samtools
* snakemake 
* whatshap (or other variant phasing tools )
### 2. Obtain the latest version of MERD 
  ```shell
     git clone https://github.com/PengJia6/MERD.git
  ```

### 3. Generate you own config file according to the config_template.yaml 

```shell
    mkdir -p /path/to/my/work/directory
    cd /path/to/my/work/directory
    cp /path/to/MERD/config_template.yaml ./config.yaml
    vim ./config.yaml # change the file according to your own requirements
```
Please see [Configration]() to check the   
### 4. Run the pipeline 
```shell
    snakemake -s /path/to/MERD/SnakeFile.smk -j -n  # check if your config.yaml is available 
    snakemake -s /path/to/MERD/SnakeFile.smk -j 4 --ri  -k # run local 
    snakemake -s /path/to/MERD/SnakeFile.smk -j 10 --ri  -k --cluster 'qsub -l nodes=1:ppn=1 -l walltime=999:00:00  # run on cluster 
```
# Configuration
```yaml
dir_work: /path/to/your/work_dir # the directory to save your output
chromsomes:
  chr1,chr2 # the chromsomes you want to process

mapping_qual_thresh: 30 # the minimum mapping quality score of read to detect recombination signal 
software: # the path of required tools 
  samtools: "/path/to/samtools" 
  whatshap: "/path/to/whatshap"
  tabix: "/path/to/tabix"

reference: /path/to/reference.genome.fa # the reference you used for reads alignments and variant calling 
family: # samples information
  family001: # unique family name 
    names: 
      paternal:
        LCL7:  # paternal sample name 
          sex: Male
      maternal:  # maternal sample name 
        LCL8
      probands: 
        - LCL5:  # name of one proband
            sex: Female
        - LCL6: name of one proband 
            sex: Female
    merged_vcf: /path/to/family1.vcf.gz # high-confidence SNVs of the family (unphased)  the samples name should be same as above 
    bams: # the bam file of the family
      LCL5:
        HiFi: /path/to/aligned_reads/family1/family1.LCL5.GRCh38.HiFi.minimap2.bam
      LCL6:
        HiFi: /path/to/aligned_reads/family1/family1.LCL6.GRCh38.HiFi.minimap2.bam
      LCL7:
        HiFi: /path/to/aligned_reads/family1/family1.LCL7.GRCh38.HiFi.minimap2.bam
      LCL8:
        HiFi: /path/to/aligned_reads/family1/family1.LCL8.GRCh38.HiFi.minimap2.bam
  HGSVC_CHS:
    names:
      paternal:
        HG00512
      maternal:
        HG00513
      probands:
        - HG00514
    merged_vcf: /path/to/HGSVC_CHS.GRCh38.HiFi.minimap2.deepvariant.raw.vcf.gz
    bams:
      HG00512:
        HiFi: /path/to/alignment/HG00512.GRCh38.HiFi.minimap2.bam
      HG00513:
        HiFi: /path/to/alignment/HG00513.GRCh38.HiFi.minimap2.bam
      HG00514:
        HiFi: /path/to/alignment/HG00514.GRCh38.HiFi.minimap2.bam
```
## Interpreting output
You may obatin a tsv file in 
**/path/to/your/work_dir/{family}{family}.{proband}.{maternal/paternal}.HiFi.final.merd.tsv** 

The columns are description as following: 
 * chrom: chromsome id of the recombination event 
 * start: start position of the recombination event
 * end: end position of the recombination event
 * Filter: The event types 
 * Primary_supports: Support reads of this event (reads with abnormal signals)
 * Average_Read_Coverage: Average read depth of the events 
 * mat_gt: genotype of maternal sample in this event
 * pat_gt: genotype of paternal sample in this event
 * proband_gt: genotype of proband sample in this event
 * total_support_re: total reads 
 * support_mut: read support the varaints phased results 
 * support_error: read does not support the varaints phased results 
 * series_support: read support not only two shift variants (this kind of reads required longer reads and were strong evidence for meiotic recombination)
 * failed_support: read does not support the shift  
 * single_support: read only support single shift (de novol mutation also had this signal)
 * shift_err: variants on this read with more than two times shifts

## Contact  
 * Please open an issue on the Github page for problems.
 * You may also contact Peng Jia at Xi'an Jiaotong University (pengjia@stu.xjtu.edu.cn)
