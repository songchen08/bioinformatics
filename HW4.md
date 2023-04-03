1. 请阐述bowtie中利用了 BWT 的什么性质提高了运算速度？并通过哪些策略优化了对内存的需求？

BWT本身并不是为mapping而开发的算法，而是一个字符串处理算法。BWT的特点在于，经过该算法的变换后，原始字符串中的相同字符会处于比较相邻的位置，从而达到压缩字符串的目的（一长串相同字符可以简化成两个字符，如将aaaaa变换成a5）。同时，BWT算法和其逆算法本质上都是逐步进行的，这使得可以对“编码”和“解码”的过程分别建立索引，通过与原字符串相同尺度的索引目录就可以快速找到某个字符在算法过程中的路径。

因此，BWT提高运算速度的性质是：由于可以跟踪字符在BWT算法中的路径，则可以快速找到原字符串中该字符上游和下游的相邻字符，这即是mapping在align序列时希望做到的事情，即BWT算法变换query序列的过程中便自然地包含了序列信息。当寻找match的位置时，便可以通过BWT的索引快速寻找小同源序列的上下游同源序列。

而BWT优化内存需求的策略则在于显著压缩了储存的基因组数据库的内存占用，而“编码”和“解码”的过程所产生的中间数据又远小于被比对序列的长度，因此节省下来的内存远远大于该算法产生的额外运算量。

2. 用bowtie将 THA2.fa mapping 到 BowtieIndex/YeastGenome 上，得到 THA2.sam，统计mapping到不同染色体上的reads数量(即统计每条染色体都map上了多少条reads)。

```
bowtie -v 2 BowtieIndex/YeastGenome -f THA2.fa -S THA2.sam
grep -v "^@" THA2.sam | cut -f 3 | sort | uniq -c
```

```
     16 chrmt
     17 chrVI
     19 chrIII
     20 chrI
     28 chrIX
     33 chrV
     50 chrII
     54 chrXI
     61 chrXIV
     69 chrVIII
     69 chrX
     72 chrXIII
     77 *
     77 chrXVI
    103 chrXV
    127 chrVII
    169 chrXII
    189 chrIV
```
3. 查阅资料，回答以下问题:

（3.1）什么是sam/bam文件中的"CIGAR string"? 它包含了什么信息?

CIGAR，Compact Idiosyncratic Gapped Alignment Report，紧凑特质间隙对齐报告，是sam/bam文件表示对齐的方式。用数字+字母来表示对齐的具体情形。M for match, X for mismatch, D for deletion, I for insertion and N for gap.

（3.2）"soft clip"的含义是什么，在CIGAR string中如何表示？

直接地说，soft clip指align不到基因组上，但是还存在于query序列比对信息之中的片段，在SEQ一列中显示（即S for soft clips）

（3.3）什么是reads的mapping quality? 它反映了什么样的信息?

为了衡量比对结果的好坏，人为定义比对质量为mapping quality，也就是MAPQ。它被定义为MAPQ=-10*log[P(mapping出错)]，即p值的延伸定义。当MAPQ较大的时候，可以认为mapping的序列是错误的，不可信的，一般将这个阈值定义为10（即有10%可能出错）。

（3.4）仅根据sam/bam文件的信息，能否推断出read mapping到的区域对应的参考基因组序列? (提示:参考https://samtools.github.io/hts-specs/SAMtags.pdf中对于MD tag的介绍)

可以。MD tag包含了所有CIGAR string中的信息，并且给出了删除、插入、非配对、gap中的具体碱基，根据这些信息可以重建出被mapping的序列的全部信息。


4. 软件安装和资源文件的下载也是生物信息学实践中的重要步骤。请自行安装教程中未涉及的bwa软件，从UCSC Genome Browser下载Yeast (S. cerevisiae, sacCer3)基因组序列。使用bwa对Yeast基因组sacCer3.fa建立索引，并利用bwa将THA2.fa，mapping到Yeast参考基因组上，并进一步转化输出得到THA2-bwa.sam文件。
```
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
tar jxf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
echo 'PATH=$PATH:/home/test/mapping/bwa-0.7.17' >> ~/.bashrc 
source ~/.bashrc #配置环境变量
./bwa  #验证是否可以成功运行

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
         pemerge       merge overlapping paired ends (EXPERIMENTAL)
         aln           gapped/ungapped alignment
         samse         generate alignment (single ended)
         sampe         generate alignment (paired ended)
         bwasw         BWA-SW for long queries

         shm           manage indices in shared memory
         fa2pac        convert FASTA to PAC format
         pac2bwt       generate BWT from PAC
         pac2bwtgen    alternative algorithm for generating BWT
         bwtupdate     update .bwt to the new format
         bwt2sa        generate SA from BWT and Occ

Note: To use BWA, you need to first index the genome with `bwa index'.
      There are three alignment algorithms in BWA: `mem', `bwasw', and
      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
      first. Please `man ./bwa.1' for the manual.

```

至此，bwa软件已经成功下载并且可以运行（需要root后使用）。从GCSC获取yeast genome文件sacCer3.fa后，
```
bwa index sacCer3.fa

# 经试验，使用mem和sw算法产生的结果与第二题不一致，后来发现这里的reads长度只有30左右，故不适合采用适用于70bp以上的mem和sw。这里考虑使用backtrack算法，使用结果第二题相差不大，以下命令均在/home/test/mapping中进行。

bwa aln sacCer3.fa THA2.fa > THA2_bwa.sai 
bwa samse sacCer3.fa THA2_bwa.sai THA2.fa > THA2_bwa.sam

# 使用第二题中的命令来表征结果并与bowtie产生结果进行对比

grep -v "^@" THA2_bwa.sam | cut -f 3 | sort | uniq -c

     24 *
     17 chrI
     54 chrII
     17 chrIII
    202 chrIV
     26 chrIX
     18 chrM
     38 chrV
     18 chrVI
    129 chrVII
     70 chrVIII
     77 chrX
     60 chrXI
    178 chrXII
     72 chrXIII
     59 chrXIV
    108 chrXV
     83 chrXVI

# 这里并没有调整参数，如果调整参数的话，或许可以获得与第二题结果更加相近甚至一样的数据。
```

extra question (in 1.1 genome browser): 利用Genome Browser浏览 1.Mapping的 Homework 得到的sam/bam文件，并仿照上文中的 examples截图展示一个 gene的区域。

图片已附在文件夹中。
